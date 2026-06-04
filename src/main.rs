use serde::Deserialize;

use std::fs::File;
use std::path::Path;
use std::env;
use std::collections::HashMap;
use std::io::BufRead;

use rayon::prelude::*;

use flate2::read::GzDecoder;

use noodles::vcf;
use noodles::vcf::header::record::value::map::{Info, Map};
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::record_buf::info::field::Value as VariantValue;
use noodles::vcf::variant::record_buf::AlternateBases as RecordAlternateBases;
use noodles::vcf::variant::record::AlternateBases; 
use noodles::vcf::variant::io::Write;

#[derive(Debug, Clone)]
pub struct StrRecord {
    pub chrom: String,
    pub pos_start: usize,
    pub pos_end: usize,
    pub motif: String,
    pub source: String,
    pub info_line: String,
}

#[derive(Debug, Deserialize, Clone)]
pub struct TandemGenotypeRow {
    #[serde(rename = "CHR")]
    pub chrom: String,
    #[serde(rename = "START")]
    pub start: usize,
    #[serde(rename = "END")]
    pub end: usize,
    #[serde(rename = "MOTIF")]
    pub motif: String,
    #[serde(rename = "GENE")]
    pub gene: String,
    #[serde(rename = "GENE_PART")]
    pub gene_part: String,
    #[serde(rename = "CN_FOR")]
    pub cn_for: String,
    #[serde(rename = "CN_REV")]
    pub cn_rev: String,
}

#[derive(Debug, Clone)]
struct TelomereCoords {
    p: (usize, usize),
    q: Option<(usize, usize)>, // Option car q peut être None (chrX, chrY)
}

#[derive(Debug, Clone)]
struct GeneInfo {
    start: usize,
    end: usize,
}

#[derive(Debug, Clone)]
struct ExonInfo {
    start: usize,
    end: usize,
}

pub fn parse_tandem_genotypes<P: AsRef<Path>>(path: P) -> Result<Vec<TandemGenotypeRow>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t') 
        .comment(Some(b'#'))
        .from_reader(file);

    let mut records = Vec::new();
    for result in rdr.deserialize() {
        let record: TandemGenotypeRow = result?;
        records.push(record);
    }
    Ok(records)
}

pub fn parse_vamos_vcf<P: AsRef<Path>>(
    path: P, 
    source_name: &str, 
    records: &mut Vec<StrRecord>
) -> Result<(), Box<dyn std::error::Error>> {
    
    let file = File::open(path)?;
    let bgzf_reader = noodles::bgzf::Reader::new(file);
    let mut reader = vcf::io::Reader::new(bgzf_reader);

    let header = reader.read_header()?; 

    for result in reader.records() {
        let record = result?;
        
        let chrom = record.reference_sequence_name().to_string();
        
        let pos_start: usize = record.variant_start()
            .ok_or("Position manquante dans le VCF")??
            .get();
        
        let motif = match record.info().get(&header, "RU") {
            Some(Ok(Some(noodles::vcf::variant::record::info::field::Value::String(ru)))) => {
                ru.split(',').next().unwrap_or(ru).to_string()
            }
            _ => {
                record.alternate_bases()
                    .iter()
                    .next()
                    .map(|alt_result| match alt_result {
                        Ok(alt) => alt.to_string(),
                        _ => "UNKNOWN".to_string()
                    })
                    .unwrap_or_else(|| "UNKNOWN".to_string())
            }
        };

        let pos_end = match record.info().get(&header, "END") {
            Some(Ok(Some(noodles::vcf::variant::record::info::field::Value::Integer(end)))) => end as usize,
            _ => pos_start + record.reference_bases().len() - 1,
        };

        records.push(StrRecord {
            chrom,
            pos_start,
            pos_end,
            motif,
            source: source_name.to_string(),
            info_line: format!("{:?}", record.info()),
        });
    }

    Ok(())
}

fn parse_genes_bed<P: AsRef<Path>>(path: P) -> Result<HashMap<String, HashMap<String, GeneInfo>>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let decoder = GzDecoder::new(file);
    let reader = std::io::BufReader::new(decoder);
    
    let mut dict_genes: HashMap<String, HashMap<String, GeneInfo>> = HashMap::new();

    for line in reader.lines() {
        let l = line?;
        if l.starts_with('#') || l.is_empty() { continue; }
        let parts: Vec<&str> = l.split('\t').collect();
        if parts.len() >= 4 {
            let chrom = parts[0].to_string();
            let start = parts[1].parse::<usize>()?;
            let end = parts[2].parse::<usize>()?;
            let gene = parts[3].to_string();

            dict_genes.entry(chrom)
                .or_insert_with(HashMap::new)
                .insert(gene, GeneInfo { start, end });
        }
    }
    Ok(dict_genes)
}

fn parse_exons_bed<P: AsRef<Path>>(path: P) -> Result<HashMap<String, HashMap<String, HashMap<String, ExonInfo>>>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let decoder = GzDecoder::new(file);
    let reader = std::io::BufReader::new(decoder);
    
    let mut dict_exons: HashMap<String, HashMap<String, HashMap<String, ExonInfo>>> = HashMap::new();

    for line in reader.lines() {
        let l = line?;
        if l.starts_with('#') || l.is_empty() { continue; }
        let parts: Vec<&str> = l.split('\t').collect();
        if parts.len() >= 4 {
            let chrom = parts[0].to_string();
            let start = parts[1].parse::<usize>()?;
            let end = parts[2].parse::<usize>()?;
            let gene_exon = parts[3];
            
            if let Some((gene, exon)) = gene_exon.split_once('_') {
                dict_exons.entry(chrom)
                    .or_insert_with(HashMap::new)
                    .entry(gene.to_string())
                    .or_insert_with(HashMap::new)
                    .insert(exon.to_string(), ExonInfo { start, end });
            }
        }
    }
    Ok(dict_exons)
}

fn build_vcf_header() -> vcf::Header {
    let mut builder = vcf::Header::builder();

    // Clé INFO pour la source
    builder = builder.add_info(
        "SOURCE",
        Map::<Info>::new(
            vcf::header::record::value::map::info::Number::Count(1),
            vcf::header::record::value::map::info::Type::String,
            "Outil d'origine ayant détecté le STR (vamos, tandem_genotypes, ou both)",
        ),
    );

    // Déclarations des métriques spécifiques à tandem-genotypes
    builder = builder.add_info(
        "TG_START",
        Map::<Info>::new(
            vcf::header::record::value::map::info::Number::Count(1),
            vcf::header::record::value::map::info::Type::Integer,
            "Position de début selon tandem-genotypes",
        ),
    );

    builder = builder.add_info(
        "TG_END",
        Map::<Info>::new(
            vcf::header::record::value::map::info::Number::Count(1),
            vcf::header::record::value::map::info::Type::Integer,
            "Position de fin selon tandem-genotypes",
        ),
    );

    builder = builder.add_info(
        "TG_MOTIF",
        Map::<Info>::new(
            vcf::header::record::value::map::info::Number::Count(1),
            vcf::header::record::value::map::info::Type::String,
            "Motif de répétition selon tandem-genotypes",
        ),
    );

    builder = builder.add_info(
        "CENTROMERE",
        Map::<Info>::new(
            vcf::header::record::value::map::info::Number::Count(0), // 0 signifie que c'est un Flag (sans valeur associée)
            vcf::header::record::value::map::info::Type::Flag,
            "Le STR intersecte une zone centromérique",
        ),
    );

    builder = builder.add_info(
        "TELOMERE_P",
        Map::<Info>::new(
            vcf::header::record::value::map::info::Number::Count(0),
            vcf::header::record::value::map::info::Type::Flag,
            "Le STR intersecte le télomère du bras p",
        ),
    );

    builder = builder.add_info(
        "TELOMERE_Q",
        Map::<Info>::new(
            vcf::header::record::value::map::info::Number::Count(0),
            vcf::header::record::value::map::info::Type::Flag,
            "Le STR intersecte le télomère du bras q",
        ),
    );
    builder = builder.add_info(
        "GENES", 
        Map::<Info>::new(
            vcf::header::record::value::map::info::Number::Count(1), 
            vcf::header::record::value::map::info::Type::String, 
            "Gènes impactés"
        )
    );
    builder = builder.add_info(
        "FEATURES", 
        Map::<Info>::new(
            vcf::header::record::value::map::info::Number::Count(1), 
            vcf::header::record::value::map::info::Type::String, 
            "Caractéristiques génomiques"
        )
    );
    builder.build()
}

fn annotate_region(
    chrom: &str, 
    start: usize, 
    end: usize, 
    centromeres: &HashMap<&str, (usize, usize)>, 
    telomeres: &HashMap<&str, TelomereCoords>,
    info_buf: &mut noodles::vcf::variant::record_buf::Info
) {
    // Vérification du Centromère
    if let Some(&(c_start, c_end)) = centromeres.get(chrom) {
        // Condition d'intersection standard entre deux intervalles [start, end] et [c_start, c_end]
        if start <= c_end && end >= c_start {
            info_buf.insert("CENTROMERE".to_string(), None); // None car c'est un Type::Flag
        }
    }

    // Vérification des Télomères
    if let Some(t_coords) = telomeres.get(chrom) {
        // Bras p
        if start <= t_coords.p.1 && end >= t_coords.p.0 {
            info_buf.insert("TELOMERE_P".to_string(), None);
        }
        // Bras q (si disponible)
        if let Some(q_coords) = t_coords.q {
            if start <= q_coords.1 && end >= q_coords.0 {
                info_buf.insert("TELOMERE_Q".to_string(), None);
            }
        }
    }
}

pub fn write_combined_vcf(
    path: &str, 
    commons: Vec<(StrRecord, TandemGenotypeRow)>, 
    vamos_only: Vec<StrRecord>, 
    tg_only: Vec<TandemGenotypeRow>,
    dict_genes: &HashMap<String, HashMap<String, GeneInfo>>,
    dict_exons: &HashMap<String, HashMap<String, HashMap<String, ExonInfo>>>
) -> Result<(), Box<dyn std::error::Error>> {
    
    let file = File::create(path)?;
    let bgzf_writer = noodles::bgzf::Writer::new(file);
    let mut writer = vcf::io::Writer::new(bgzf_writer);
    
    let header = build_vcf_header();
    writer.write_header(&header)?;

    let centromere_coords: HashMap<&str, (usize, usize)> = HashMap::from([
        ("chr1", (121535434, 124535434)), ("chr2", (92326171, 95326171)),
        ("chr3", (90504854, 93504854)),   ("chr4", (49660117, 52660117)),
        ("chr5", (46405641, 49405641)),   ("chr6", (58626368, 61626368)),
        ("chr7", (58169654, 60828234)),   ("chr8", (44033745, 45877265)),
        ("chr9", (43236168, 45518558)),   ("chr10", (39686683, 41593521)),
        ("chr11", (51078349, 54425074)), ("chr12", (34769408, 37185252)),
        ("chr13", (16000001, 18051248)), ("chr14", (16000001, 18173523)),
        ("chr15", (17000001, 19725254)), ("chr16", (36311159, 38280682)),
        ("chr17", (22813680, 26885980)), ("chr18", (15460900, 20861206)),
        ("chr19", (24498981, 27190874)), ("chr20", (26436233, 30038348)),
        ("chr21", (10864561, 12915808)), ("chr22", (12954789, 15054318)),
        ("chrX", (58605580, 62412542)),   ("chrY", (10316945, 10544039)),
    ]);

    let hg38_telomeres: HashMap<&str, TelomereCoords> = HashMap::from([
        ("chr1",  TelomereCoords { p: (0, 10000), q: Some((248946422, 248956422)) }),
        ("chr2",  TelomereCoords { p: (0, 10000), q: Some((242183529, 242193529)) }),
        ("chr3",  TelomereCoords { p: (0, 10000), q: Some((198285559, 198295559)) }),
        ("chr4",  TelomereCoords { p: (0, 10000), q: Some((190204555, 190214555)) }),
        ("chr5",  TelomereCoords { p: (0, 10000), q: Some((181528259, 181538259)) }),
        ("chr6",  TelomereCoords { p: (0, 10000), q: Some((170795979, 170805979)) }),
        ("chr7",  TelomereCoords { p: (0, 10000), q: Some((159335973, 159345973)) }),
        ("chr8",  TelomereCoords { p: (0, 10000), q: Some((145128636, 145138636)) }),
        ("chr9",  TelomereCoords { p: (0, 10000), q: Some((138384717, 138394717)) }),
        ("chr10", TelomereCoords { p: (0, 10000), q: Some((133787422, 133797422)) }),
        ("chr11", TelomereCoords { p: (0, 10000), q: Some((135076622, 135086622)) }),
        ("chr12", TelomereCoords { p: (0, 10000), q: Some((133265309, 133275309)) }),
        ("chr13", TelomereCoords { p: (0, 10000), q: Some((114354328, 114364328)) }),
        ("chr14", TelomereCoords { p: (0, 10000), q: Some((107033718, 107043718)) }),
        ("chr15", TelomereCoords { p: (0, 10000), q: Some((101981189, 101991189)) }),
        ("chr16", TelomereCoords { p: (0, 10000), q: Some((90328345,  90338345))  }),
        ("chr17", TelomereCoords { p: (0, 10000), q: Some((83247441,  83257441))  }),
        ("chr18", TelomereCoords { p: (0, 10000), q: Some((80363285,  80373285))  }),
        ("chr19", TelomereCoords { p: (0, 10000), q: Some((58607616,  58617616))  }),
        ("chr20", TelomereCoords { p: (0, 10000), q: Some((64434167,  64444167))  }),
        ("chr21", TelomereCoords { p: (0, 10000), q: Some((46699983,  46709983))  }),
        ("chr22", TelomereCoords { p: (0, 10000), q: Some((50808468,  50818468))  }),
        ("chrX",  TelomereCoords { p: (0, 10000), q: None }),
        ("chrY",  TelomereCoords { p: (0, 10000), q: None }),
    ]);

    // ---- 1. PARALLÉLISATION : Génération des Records VCF en mémoire ----
    
    // On transforme le vecteur `commons` en parallèle grâce à .into_par_iter()
    let records_commons: Vec<RecordBuf> = commons.into_par_iter().map(|(v_rec, tg_rec)| {
        let vcf_pos = noodles::core::Position::try_from(v_rec.pos_start).unwrap();

        let mut record_builder = RecordBuf::builder()
            .set_reference_sequence_name(v_rec.chrom.clone())
            .set_variant_start(vcf_pos) 
            .set_reference_bases("N") 
            .set_alternate_bases(RecordAlternateBases::from(vec![String::from("<STR>")]));

        let mut info_buf = noodles::vcf::variant::record_buf::Info::default();
        info_buf.insert("SOURCE".to_string(), Some(VariantValue::from("both")));
        info_buf.insert("TG_START".to_string(), Some(VariantValue::from(tg_rec.start as i32)));
        info_buf.insert("TG_END".to_string(), Some(VariantValue::from(tg_rec.end as i32)));
        info_buf.insert("TG_MOTIF".to_string(), Some(VariantValue::from(tg_rec.motif.clone())));
        info_buf.insert("RU".to_string(), Some(VariantValue::from(v_rec.motif.clone())));

        let final_start = std::cmp::min(v_rec.pos_start, tg_rec.start);
        let final_end = std::cmp::max(v_rec.pos_end, tg_rec.end);
        
        annotate_regions_full(&v_rec.chrom, final_start, final_end, dict_genes, dict_exons, &centromere_coords, &hg38_telomeres, &mut info_buf);
        
        record_builder.set_info(info_buf).build()
    }).collect();

    // On fait de même pour `vamos_only`
    let records_vamos: Vec<RecordBuf> = vamos_only.into_par_iter().map(|v_rec| {
        let vcf_pos = noodles::core::Position::try_from(v_rec.pos_start).unwrap();

        let mut info_buf = noodles::vcf::variant::record_buf::Info::default();
        info_buf.insert("SOURCE".to_string(), Some(VariantValue::from(v_rec.source.clone())));
        info_buf.insert("RU".to_string(), Some(VariantValue::from(v_rec.motif.clone())));

        annotate_regions_full(&v_rec.chrom, v_rec.pos_start, v_rec.pos_end, dict_genes, dict_exons, &centromere_coords, &hg38_telomeres, &mut info_buf);

        RecordBuf::builder()
            .set_reference_sequence_name(v_rec.chrom)
            .set_variant_start(vcf_pos)
            .set_reference_bases("N")
            .set_alternate_bases(RecordAlternateBases::from(vec![String::from("<STR>")]))
            .set_info(info_buf)
            .build()
    }).collect();

    // On fait de même pour `tg_only`
    let records_tg: Vec<RecordBuf> = tg_only.into_par_iter().map(|tg_rec| {
        let vcf_pos = noodles::core::Position::try_from(tg_rec.start).unwrap(); 

        let mut info_buf = noodles::vcf::variant::record_buf::Info::default();
        info_buf.insert("SOURCE".to_string(), Some(VariantValue::from("tandem_genotypes")));
        info_buf.insert("TG_START".to_string(), Some(VariantValue::from(tg_rec.start as i32)));
        info_buf.insert("TG_END".to_string(), Some(VariantValue::from(tg_rec.end as i32)));
        info_buf.insert("TG_MOTIF".to_string(), Some(VariantValue::from(tg_rec.motif.clone())));

        annotate_regions_full(&tg_rec.chrom, tg_rec.start, tg_rec.end, dict_genes, dict_exons, &centromere_coords, &hg38_telomeres, &mut info_buf);

        RecordBuf::builder()
            .set_reference_sequence_name(tg_rec.chrom)
            .set_variant_start(vcf_pos) 
            .set_reference_bases("N")
            .set_alternate_bases(RecordAlternateBases::from(vec![String::from("<STR>")]))
            .set_info(info_buf)
            .build()
    }).collect();


    // ---- 2. ÉCRITURE SÉQUENTIELLE DANS LE FICHIER VCF ----
    // Rayon garantit que l'ordre initial des éléments est préservé lors du .collect()
    
    for record in records_commons {
        writer.write_variant_record(&header, &record)?;
    }
    for record in records_vamos {
        writer.write_variant_record(&header, &record)?;
    }
    for record in records_tg {
        writer.write_variant_record(&header, &record)?;
    }

    let bgzf_writer = writer.into_inner();
    bgzf_writer.finish()?;

    Ok(())
}

fn annotate_regions_full(
    chrom: &str,
    pos: usize,
    pos_end: usize,
    dict_genes: &HashMap<String, HashMap<String, GeneInfo>>,
    dict_exons: &HashMap<String, HashMap<String, HashMap<String, ExonInfo>>>,
    centromeres: &HashMap<&str, (usize, usize)>,
    telomeres: &HashMap<&str, TelomereCoords>,
    info_buf: &mut noodles::vcf::variant::record_buf::Info,
) {
    let mut list_genes = Vec::new();
    let mut list_features = Vec::new();

    // 1. Analyse des gènes et exons
    if let Some(genes_on_chrom) = dict_genes.get(chrom) {
        for (g, gene_info) in genes_on_chrom {
            if pos_end >= gene_info.start && pos <= gene_info.end {
                list_genes.push(g.clone());

                if let Some(chrom_exons) = dict_exons.get(chrom) {
                    if let Some(exons) = chrom_exons.get(g) {
                        if !exons.is_empty() {
                            // Trier les exons numériquement (ex: "exon1", "exon2")
                            let mut sorted_exons: Vec<&String> = exons.keys().collect();
                            sorted_exons.sort_by_key(|key| {
                                key.replace("exon", "").parse::<usize>().unwrap_or(0)
                            });

                            let first_exon_name = sorted_exons[0];
                            let last_exon_name = sorted_exons[sorted_exons.len() - 1];

                            let first_exon_start = exons[first_exon_name].start;
                            let last_exon_start = exons[last_exon_name].start;

                            let strand = if first_exon_start < last_exon_start { "+" } else { "-" };

                            // Flag pour savoir si on a déclenché un UTR
                            let mut utr_triggered = false;

                            // Vérification 5' UTR
                            if (strand == "+" && pos < exons[first_exon_name].start)
                                || (strand == "-" && pos > exons[first_exon_name].end)
                            {
                                list_features.push("5'UTR".to_string());
                                utr_triggered = true;
                            }
                            // Vérification 3' UTR
                            else if (strand == "+" && pos > exons[last_exon_name].end)
                                || (strand == "-" && pos < exons[last_exon_name].start)
                            {
                                list_features.push("3'UTR".to_string());
                                utr_triggered = true;
                            }

                            if !utr_triggered {
                                let mut found_exon = false;
                                // Parcourir tous les exons pour voir s'il est dedans
                                for (e, exon_info) in exons {
                                    if exon_info.start < pos && pos < exon_info.end {
                                        list_features.push(e.clone());
                                        found_exon = true;
                                    }
                                }

                                // Si ce n'est ni UTR ni exonic, c'est intronique
                                if !found_exon {
                                    list_features.push("intronic".to_string());
                                }
                            }
                        } else {
                            list_features.push("intronic".to_string());
                        }
                    } else {
                        list_features.push("intronic".to_string());
                    }
                } else {
                    list_features.push("intronic".to_string());
                }
            }
        }
    }

    // 2. Centromères
    if let Some(&(c_start, c_end)) = centromeres.get(chrom) {
        if pos > c_start && pos < c_end {
            list_features.push("centromere".to_string());
        }
    }

    // 3. Télomères
    if let Some(t_coords) = telomeres.get(chrom) {
        if let Some(q_coords) = t_coords.q {
            if pos > q_coords.0 && pos < q_coords.1 {
                list_features.push("telomere_q".to_string());
            }
        }
        if pos > t_coords.p.0 && pos < t_coords.p.1 {
            list_features.push("telomere_p".to_string());
        }
    }

    // 4. Traitement des cas vides par défaut
    if list_genes.is_empty() {
        list_genes.push("intergenic".to_string());
    }
    if list_features.is_empty() {
        list_features.push(".".to_string());
    }

    // Enlever les doublons potentiels dans les features
    list_features.dedup();

    // 5. Injection dans le buffer INFO du VCF Noodles
    info_buf.insert("GENES".to_string(), Some(VariantValue::from(list_genes.join(","))));
    info_buf.insert("FEATURES".to_string(), Some(VariantValue::from(list_features.join(","))));
}

fn normalize_motif(motif: &str) -> String {
    let motif = motif.to_uppercase();
    let len = motif.len();
    let mut min_rotation = motif.clone();
    
    for i in 0..len {
        let rotation = format!("{}{}", &motif[i..], &motif[0..i]);
        if rotation < min_rotation {
            min_rotation = rotation;
        }
    }
    min_rotation
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let margin = 50; 
    let args: Vec<String> = env::args().collect();
    
    let dict_genes = parse_genes_bed("/home/brandaoq/scripts/genes.bed.gz")?;
    let dict_exons = parse_exons_bed("/home/brandaoq/scripts/MANE_Select_exons.bed.gz")?;
    // Vérification de sécurité pour les arguments
    if args.len() < 5 {
        eprintln!("Usage: {} <vamos_hap1.vcf> <vamos_hap2.vcf> <tg_records.tsv> <output_file>", args[0]);
        std::process::exit(1);
    }

    // 1. Lire VaMoS Hap1 et Hap2 dans une liste temporaire
    let mut raw_vamos_records = Vec::new();
    parse_vamos_vcf(&args[1], "vamos_hap1", &mut raw_vamos_records)?;
    parse_vamos_vcf(&args[2], "vamos_hap2", &mut raw_vamos_records)?;
    
    println!("Nombre brut de STRs lus chez VaMoS (Hap1 + Hap2) : {}", raw_vamos_records.len());

    // ---- FUSION DES DOUBLONS (Hap1 & Hap2) EN AMONT ----
    let mut collapsed_vamos: HashMap<(String, usize, String), StrRecord> = HashMap::new();

    for mut record in raw_vamos_records {
        // On crée une clé unique basée sur : Chromosome, Position de début, et le Motif Normalisé
        let key = (
            record.chrom.clone(),
            record.pos_start,
            normalize_motif(&record.motif),
        );

        // Si la clé existe déjà, c'est que le variant est présent dans Hap1 ET Hap2
        if let Some(existing_record) = collapsed_vamos.get_mut(&key) {
            existing_record.source = "vamos_hap1_and_hap2".to_string();
        } else {
            // Sinon, on l'ajoute pour la première fois
            collapsed_vamos.insert(key, record);
        }
    }

    // On transforme notre dictionnaire fusionné en un vecteur pour la suite de l'algorithme
    let all_vamos_records: Vec<StrRecord> = collapsed_vamos.into_values().collect();
    println!("Nombre de STRs VaMoS uniques après fusion des haplotypes : {}", all_vamos_records.len());

    // 2. Lire Tandem-Genotypes TSV
    let tg_records = parse_tandem_genotypes(&args[3])?;

    // 3. Classification
    let mut commons = Vec::new();
    let mut vamos_only = Vec::new();
    let mut tg_only = tg_records.clone(); 

    // On itère sur la liste dédoublonnée
    for v_rec in all_vamos_records {
        let mut found = false;
        
        for idx in 0..tg_only.len() {
            let tg_rec = &tg_only[idx];
            
            // Calcul de la distance absolue de manière sécurisée en Rust
            let diff_start = (v_rec.pos_start as isize - tg_rec.start as isize).abs();
            
            if v_rec.chrom == tg_rec.chrom 
               && diff_start <= margin
               && normalize_motif(&v_rec.motif) == normalize_motif(&tg_rec.motif) 
            {
                commons.push((v_rec.clone(), tg_rec.clone()));
                tg_only.remove(idx); // Retire l'élément trouvé
                found = true;
                break;
            }
        }
        
        if !found {
            vamos_only.push(v_rec);
        }
    }

    println!("STR Communs (VaMoS + Tandem-Genotypes) : {}", commons.len());
    println!("STR VaMoS uniques : {}", vamos_only.len());
    println!("STR Tandem-Genotypes uniques : {}", tg_only.len());

    write_combined_vcf(&args[4], commons, vamos_only, tg_only, &dict_genes, &dict_exons)?;

    Ok(())
}