use serde::Deserialize;

use std::fs::File;
use std::path::Path;
use std::env;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Seek}; // Ajout de Read et Seek pour la détection magique

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
    pub chrom: String,     // Colonne 1 (index 0) : ex: chr1
    pub start: usize,      // Colonne 2 (index 1) : ex: 231799889
    pub end: usize,        // Colonne 3 (index 2) : ex: 231799934
    pub motif: String,     // Colonne 4 (index 3) : ex: TCCCTTCCTCCCTTCC
    pub gene: String,      // Colonne 5 (index 4) : Nom du gène
    pub gene_part: String, // Colonne 6 (index 5) : Exon/Intron...
    pub cn_for: String,    // Colonne 7 (index 6) : Copy number forward
    pub cn_rev: String,    // Colonne 8 (index 7) : Copy number reverse
}

#[derive(Debug, Clone)]
struct TelomereCoords {
    p: (usize, usize),
    q: Option<(usize, usize)>, 
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

/// Fonction utilitaire ouvrant un fichier de manière transparente, 
/// qu'il soit compressé (Gzip/BGZF) ou en texte brut.
fn open_vcf_file<P: AsRef<Path>>(path: P) -> Result<Box<dyn BufRead>, Box<dyn std::error::Error>> {
    let mut file = File::open(path)?;
    
    // On inspecte les 2 premiers octets (Magic Number)
    let mut header = [0u8; 2];
    let bytes_read = file.read(&mut header)?;
    
    // On rembobine pour que Noodles lise dès le début
    file.rewind()?;

    if bytes_read == 2 && header == [0x1f, 0x8b] {
        // Format compressé
        let decoder = GzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        // Format texte brut
        Ok(Box::new(BufReader::new(file)))
    }
}

pub fn parse_tandem_genotypes<P: AsRef<Path>>(path: P) -> Result<Vec<TandemGenotypeRow>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let mut records = Vec::new();

    // On crée le parser CSV configuré sans en-tête et tolérant
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false) // Traite la première ligne directement comme de la donnée
        .flexible(true)     // Évite les erreurs si le nombre de colonnes varie
        .from_reader(std::io::Cursor::new(Vec::new())); // Instance temporaire pour la configuration

    for line_result in reader.lines() {
        let line = line_result?;
        // Ignore manuellement toutes les lignes de commentaires ou les lignes vides
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        // On parse la ligne brute tabulée de manière sécurisée
        let mut string_reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .flexible(true)
            .from_reader(line.as_bytes());

        if let Some(result) = string_reader.deserialize::<TandemGenotypeRow>().next() {
            match result {
                Ok(record) => records.push(record),
                Err(e) => {
                    eprintln!("Avertissement : Impossible de lire une ligne Tandem-Genotypes : {}", e);
                }
            }
        }
    }

    Ok(records)
}

pub fn parse_vamos_vcf<P: AsRef<Path>>(
    path: P, 
    source_name: &str, 
    records: &mut Vec<StrRecord>
) -> Result<(), Box<dyn std::error::Error>> {
    
    // Utilisation du lecteur universel (gère texte brut ET g釐/bgzf)
    let universal_reader = open_vcf_file(path)?;
    let mut reader = vcf::io::Reader::new(universal_reader);

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

    builder = builder.add_info(
        "SOURCE",
        Map::<Info>::new(
            vcf::header::record::value::map::info::Number::Count(1),
            vcf::header::record::value::map::info::Type::String,
            "Outil d'origine ayant détecté le STR (vamos, tandem_genotypes, ou both)",
        ),
    );

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
            vcf::header::record::value::map::info::Number::Count(0), 
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
    builder = builder.add_info(
        "CN_FOR",
        Map::<Info>::new(
            vcf::header::record::value::map::info::Number::Unknown,
            vcf::header::record::value::map::info::Type::String,
            "Copy number change in each DNA read covering the forward strand",
        ),
    );
    builder = builder.add_info(
        "CN_REV",
        Map::<Info>::new(
            vcf::header::record::value::map::info::Number::Unknown,
            vcf::header::record::value::map::info::Type::String,
            "Copy number change in each DNA read covering the reverse strand",
        ),
    );
    builder.build()
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
        info_buf.insert("CN_FOR".to_string(), Some(VariantValue::from(tg_rec.cn_for.clone())));
        info_buf.insert("CN_REV".to_string(), Some(VariantValue::from(tg_rec.cn_rev.clone())));

        let final_start = std::cmp::min(v_rec.pos_start, tg_rec.start);
        let final_end = std::cmp::max(v_rec.pos_end, tg_rec.end);
        
        annotate_regions_full(&v_rec.chrom, final_start, final_end, dict_genes, dict_exons, &centromere_coords, &hg38_telomeres, &mut info_buf);
        
        record_builder.set_info(info_buf).build()
    }).collect();

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

    let records_tg: Vec<RecordBuf> = tg_only.into_par_iter().map(|tg_rec| {
        let vcf_pos = noodles::core::Position::try_from(tg_rec.start).unwrap(); 

        let mut info_buf = noodles::vcf::variant::record_buf::Info::default();
        info_buf.insert("SOURCE".to_string(), Some(VariantValue::from("tandem_genotypes")));
        info_buf.insert("TG_START".to_string(), Some(VariantValue::from(tg_rec.start as i32)));
        info_buf.insert("TG_END".to_string(), Some(VariantValue::from(tg_rec.end as i32)));
        info_buf.insert("TG_MOTIF".to_string(), Some(VariantValue::from(tg_rec.motif.clone())));
        info_buf.insert("CN_FOR".to_string(), Some(VariantValue::from(tg_rec.cn_for.clone())));
        info_buf.insert("CN_REV".to_string(), Some(VariantValue::from(tg_rec.cn_rev.clone())));
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
            // Utilisation robuste de l'intervalle [pos, pos_end]
            if pos_end >= gene_info.start && pos <= gene_info.end {
                list_genes.push(g.clone());

                if let Some(chrom_exons) = dict_exons.get(chrom) {
                    if let Some(exons) = chrom_exons.get(g) {
                        if !exons.is_empty() {
                            let mut sorted_exons: Vec<&String> = exons.keys().collect();
                            sorted_exons.sort_by_key(|key| {
                                key.replace("exon", "").parse::<usize>().unwrap_or(0)
                            });

                            let first_exon_name = sorted_exons[0];
                            let last_exon_name = sorted_exons[sorted_exons.len() - 1];

                            let first_exon_start = exons[first_exon_name].start;
                            let last_exon_start = exons[last_exon_name].start;

                            let strand = if first_exon_start < last_exon_start { "+" } else { "-" };
                            let mut utr_triggered = false;

                            // Vérification 5' UTR avec prise en compte de la largeur du variant
                            if (strand == "+" && pos < exons[first_exon_name].start)
                                || (strand == "-" && pos_end > exons[first_exon_name].end)
                            {
                                list_features.push("5'UTR".to_string());
                                utr_triggered = true;
                            }
                            // Vérification 3' UTR avec prise en compte de la largeur du variant
                            else if (strand == "+" && pos_end > exons[last_exon_name].end)
                                || (strand == "-" && pos < exons[last_exon_name].start)
                            {
                                list_features.push("3'UTR".to_string());
                                utr_triggered = true;
                            }

                            if !utr_triggered {
                                let mut found_exon = false;
                                for (e, exon_info) in exons {
                                    // Intersection d'intervalles propre
                                    if pos_end >= exon_info.start && pos <= exon_info.end {
                                        list_features.push(e.clone());
                                        found_exon = true;
                                    }
                                }

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

    // 2. Centromères avec intersection propre
    if let Some(&(c_start, c_end)) = centromeres.get(chrom) {
        if pos_end >= c_start && pos <= c_end {
            list_features.push("centromere".to_string());
        }
    }

    // 3. Télomères avec intersection propre
    if let Some(t_coords) = telomeres.get(chrom) {
        if let Some(q_coords) = t_coords.q {
            if pos_end >= q_coords.0 && pos <= q_coords.1 {
                list_features.push("telomere_q".to_string());
            }
        }
        if pos_end >= t_coords.p.0 && pos <= t_coords.p.1 {
            list_features.push("telomere_p".to_string());
        }
    }

    if list_genes.is_empty() {
        list_genes.push("intergenic".to_string());
    }
    if list_features.is_empty() {
        list_features.push(".".to_string());
    }

    list_features.dedup();

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
    
    if args.len() < 5 {
        eprintln!("Usage: {} <vamos_hap1.vcf> <vamos_hap2.vcf> <tg_records.tsv> <output_file>", args[0]);
        std::process::exit(1);
    }

    let dict_genes = parse_genes_bed("/home/brandaoq/scripts/genes.bed.gz")?;
    let dict_exons = parse_exons_bed("/home/brandaoq/scripts/MANE_Select_exons.bed.gz")?;

    let mut raw_vamos_records = Vec::new();
    parse_vamos_vcf(&args[1], "vamos_hap1", &mut raw_vamos_records)?;
    parse_vamos_vcf(&args[2], "vamos_hap2", &mut raw_vamos_records)?;
    
    println!("Nombre brut de STRs lus chez VaMoS (Hap1 + Hap2) : {}", raw_vamos_records.len());

    let mut collapsed_vamos: HashMap<(String, usize, String), StrRecord> = HashMap::new();

    // Suppression du 'mut' inutile devant 'record' pour éliminer l'avertissement
    for record in raw_vamos_records {
        let key = (
            record.chrom.clone(),
            record.pos_start,
            normalize_motif(&record.motif),
        );

        if let Some(existing_record) = collapsed_vamos.get_mut(&key) {
            existing_record.source = "vamos_hap1_and_hap2".to_string();
        } else {
            collapsed_vamos.insert(key, record);
        }
    }

    let all_vamos_records: Vec<StrRecord> = collapsed_vamos.into_values().collect();
    println!("Nombre de STRs VaMoS uniques après fusion des haplotypes : {}", all_vamos_records.len());

    let tg_records = parse_tandem_genotypes(&args[3])?;

    let mut commons = Vec::new();
    let mut vamos_only = Vec::new();
    let mut tg_only = tg_records.clone(); 

    for v_rec in all_vamos_records {
        let mut found = false;
        
        for idx in 0..tg_only.len() {
            let tg_rec = &tg_only[idx];
            let diff_start = (v_rec.pos_start as isize - tg_rec.start as isize).abs();
            
            if v_rec.chrom == tg_rec.chrom 
               && diff_start <= margin
               && normalize_motif(&v_rec.motif) == normalize_motif(&tg_rec.motif) 
            {
                commons.push((v_rec.clone(), tg_rec.clone()));
                tg_only.remove(idx); 
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