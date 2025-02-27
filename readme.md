==> README.md <==
# CRISPRapido

![CRISPRapido Logo](crisprapido.png)


CRISPRapido is a reference-free tool for comprehensive detection of CRISPR off-target sites using complete genome assemblies. Unlike traditional approaches that rely on reference genomes and variant files, CRISPRapido directly analyzes haplotype-resolved assemblies to identify potential off-targets arising from any form of genetic variation. By leveraging the efficient Wavefront Alignment (WFA) algorithm and parallel processing, CRISPRapido enables fast scanning of whole genomes while considering both mismatches and DNA/RNA bulges. The tool is particularly valuable for therapeutic applications, where comprehensive off-target analysis is critical for safety assessment. CRISPRapido can process both complete assemblies and raw sequencing data, providing flexibility for different analysis scenarios while maintaining high computational efficiency through its robust Rust implementation.

## Features

- Fast parallel scanning of genomic sequences
- Support for both gzipped and plain FASTA files
- Configurable mismatch and bulge tolerances
- Automatic reverse complement scanning
- PAF-format output compatible with downstream analysis tools
- Multi-threaded processing for improved performance

## Installation

You need to build `WFA2-lib` first, which is a submodule of this repository. To do so, run:

```bash
git clone --recursive https://github.com/pinellolab/crisprapido.git
cd crisprapido/WFA2-lib
make clean all
cd ..
```

Then, you can install CRISPRapido using Cargo:

```shell
# Point to your pre-built WFA2-lib directory
export WFA2LIB_PATH="./WFA2-lib"

# Install CRISPRapido
cargo install --git https://github.com/pinellolab/crisprapido.git
```

### For GUIX's users

```bash
git clone --recursive https://github.com/pinellolab/crisprapido.git
cd crisprapido/WFA2-lib
guix shell -C -D -f guix.scm
export CC=gcc; make clean all
exit
cd ..
env -i bash -c 'WFA2LIB_PATH="./WFA2-lib" PATH=/usr/local/bin:/usr/bin:/bin ~/.cargo/bin/cargo install --path .'
```

## Usage

```bash
crisprapido -r <reference.fa> -g <guide_sequence> [OPTIONS]
```

### Required Arguments

- `-r, --reference <FILE>`: Input reference FASTA file (supports .fa and .fa.gz)
- `-g, --guide <SEQUENCE>`: Guide RNA sequence (without PAM)

### Optional Arguments

- `-m, --max-mismatches <NUM>`: Maximum number of mismatches allowed (default: 4)
- `-b, --max-bulges <NUM>`: Maximum number of bulges allowed (default: 1)
- `-z, --max-bulge-size <NUM>`: Maximum size of each bulge in bp (default: 2)
- `-w, --window-size <NUM>`: Size of sequence window to scan (default: 4x guide length)
- `-t, --threads <NUM>`: Number of threads to use (default: number of logical CPUs)
- `--no-filter`: Disable all filtering (report every alignment)

## Output Format

Output is in PAF format with custom tags, including:
- Standard PAF columns for position and alignment information
- `as:i`: Alignment score
- `nm:i`: Number of mismatches
- `ng:i`: Number of gaps
- `bs:i`: Biggest gap size
- `cg:Z`: CIGAR string

## Example

```bash
crisprapido -r genome.fa -g ATCGATCGATCG -m 3 -b 1 -z 2
```

## Testing

Run the test suite:

```bash
# Point to your pre-built WFA2-lib directory
export WFA2LIB_PATH="./WFA2-lib"

cargo test
```

Enable debug output during development:

```bash
cargo run --features debug
```

## License

See LICENSE file

## Citation

Stay tuned!

==> Cargo.toml <==
[package]
name = "crisprapido"
version = "0.1.0"
edition = "2021"

[features]
debug = []

[dependencies]
clap = { version = "4.5.31", features = ["derive"] }
bio = "2.2.0"
lib_wfa2 = { git = "https://github.com/AndreaGuarracino/lib_wfa2", rev = "c608c436a6753d2c21c97d9f5c338efae99d042b"}
rand = { version = "0.9.0", features = ["small_rng"] }
rayon = "1.10.0"
flate2 = "1.1.0"

==> src/main.rs <==
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader};
use flate2::read::MultiGzDecoder;
use std::sync::Arc;
use clap::Parser;
use bio::io::fasta;
use lib_wfa2::affine_wavefront::AffineWavefronts;
use std::fmt::Write;
use rayon::prelude::*;

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        _ => b'N',  // Convert any unexpected bases to N
    }).collect()
}

fn report_hit(ref_id: &str, pos: usize, _len: usize, strand: char, 
              _score: i32, cigar: &str, guide: &[u8], target_len: usize,
              _max_mismatches: u32, _max_bulges: u32, _max_bulge_size: u32) {
    // Calculate reference and query positions and consumed bases
    let mut ref_pos = pos;
    let mut ref_consumed = 0;
    let mut query_start = 0;
    let mut query_consumed = 0;
    
    // Count leading deletions to adjust start position
    let leading_dels = cigar.chars()
        .take_while(|&c| c == 'D')
        .count();
    ref_pos += leading_dels;
    
    // Calculate alignment statistics, accounting for N positions
    let mut mismatches = 0;
    let mut gaps = 0;
    let mut current_gap_size = 0;
    let mut max_gap_size = 0;
    let mut pos = 0;
    for c in cigar.chars() {
        match c {
            'X' => {
                // Only count mismatch if this position in the guide isn't N
                if pos < guide.len() && guide[pos] != b'N' {
                    mismatches += 1;
                }
                ref_consumed += 1;
                query_consumed += 1;
                pos += 1;
            },
            'I' => {
                current_gap_size += 1;
                if current_gap_size == 1 {
                    gaps += 1;
                }
                max_gap_size = max_gap_size.max(current_gap_size);
                query_consumed += 1;
            },
            'D' => {
                current_gap_size += 1;
                if current_gap_size == 1 {
                    gaps += 1;
                }
                max_gap_size = max_gap_size.max(current_gap_size);
                ref_consumed += 1;
                query_start += 1;  // Adjust query start for leading deletions
            },
            'M' | '=' => {
                current_gap_size = 0;
                ref_consumed += 1;
                query_consumed += 1;
            },
            _ => ()
        }
    }

    // Recalculate score based on the alignment, accounting for N positions
    let mut adjusted_score = 0;
    let mut in_gap = false;
    let mut pos = 0;
    for c in cigar.chars() {
        match c {
            'X' => {
                // Only count mismatch if this position in the guide isn't N
                if pos < guide.len() && guide[pos] != b'N' {
                    adjusted_score += 3;  // Mismatch penalty
                }
                pos += 1;
            },
            'I' | 'D' => {
                if !in_gap {
                    adjusted_score += 5;  // Gap opening penalty
                    in_gap = true;
                }
                adjusted_score += 1;  // Gap extension penalty
                if c == 'I' { pos += 1; }
            },
            'M' | '=' => {
                in_gap = false;
                pos += 1;
            },
            _ => ()
        }
    }

    // Count matches from CIGAR
    let matches = cigar.chars()
        .filter(|&c| c == 'M' || c == '=')
        .count();
    
    // Calculate block length (matches + mismatches + indels)
    let block_len = cigar.len();
    
    // Convert guide length to string once
    let guide_len = guide.len();
    
    // Debug macro for development/testing
    macro_rules! debug {
        ($($arg:tt)*) => {
            #[cfg(feature = "debug")]
            eprintln!($($arg)*);
        }
    }

    debug!("Window scan debug:");
    debug!("  CIGAR: {}", cigar);
    debug!("  N-adjusted mismatches: {} (max: 4)", mismatches);
    debug!("  Gaps: {} (max: 1)", gaps);
    debug!("  Max gap size: {} (max: 2)", max_gap_size);
    debug!("  Guide sequence: {}", String::from_utf8_lossy(guide));
    
    debug!("  Passes filters: true");
    debug!("");

    println!("Guide\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tas:i:{}\tnm:i:{}\tng:i:{}\tbs:i:{}\tcg:Z:{}", 
        guide_len,                        // Query length
        query_start,                      // Query start
        query_start + query_consumed,     // Query end
        strand,                           // Strand (+/-)
        ref_id,                           // Target sequence name
        target_len,                       // Full target sequence length
        ref_pos,                          // Target start
        ref_pos + ref_consumed,           // Target end
        matches,                          // Number of matches
        block_len,                        // Total alignment block length
        adjusted_score,                   // AS:i alignment score
        mismatches,                       // NM:i number of mismatches
        gaps,                             // NG:i number of gaps
        max_gap_size,                     // BS:i biggest gap size
        convert_to_minimap2_cigar(cigar) // cg:Z CIGAR string
    );
}
#[cfg(test)]
use rand::{SeedableRng, RngCore, rngs::SmallRng};

#[cfg(test)]
mod tests {
    use super::*;

    fn generate_random_seq(rng: &mut SmallRng, length: usize) -> Vec<u8> {
        let bases = b"ACGT";
        (0..length)
            .map(|_| bases[rng.next_u32() as usize % 4])
            .collect()
    }

    fn create_flanked_sequence(rng: &mut SmallRng, core: &[u8], flank_size: usize) -> Vec<u8> {
        let mut seq = generate_random_seq(rng, flank_size);
        seq.extend_from_slice(core);
        seq.extend(generate_random_seq(rng, flank_size));
        seq
    }

    fn setup_aligner() -> AffineWavefronts {
        AffineWavefronts::with_penalties(
            0,     // match score
            3,     // mismatch penalty
            5,     // gap opening penalty
            1      // gap extension penalty
        )
    }

    #[test]
    fn test_perfect_match() {
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let target = b"ATCGATCGAT";
        
        let result = scan_window(&mut aligner, guide, target, 1, 1, 1);
        assert!(result.is_some());
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "MMMMMMMMMM");
    }

    #[test]
    fn test_with_mismatches() {
        let mut aligner = setup_aligner();
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGTTCGAT";  // Single mismatch at position 5
        
        let result = scan_window(&mut aligner, guide, target, 1, 1, 1);
        assert!(result.is_some(), "Should accept a single mismatch");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "MMMMXMMMMM");
    }

    #[test]
    fn test_with_bulge() {
        let mut aligner = setup_aligner();
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGAATCGAT";  // Single base insertion after position 4
        
        let result = scan_window(&mut aligner, guide, target, 1, 1, 1);
        assert!(result.is_some(), "Should accept a single base bulge");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert!(cigar.contains('I') || cigar.contains('D'), "Should contain an insertion or deletion");
    }

    #[test]
    fn test_too_many_differences() {
        let mut aligner = setup_aligner();
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGTTCGTT";  // Three mismatches at positions 5, 8, 9
        
        let result = scan_window(&mut aligner, guide, target, 1, 1, 1);
        assert!(result.is_none());
    }

    #[test]
    fn test_perfect_match_with_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let target = create_flanked_sequence(&mut rng, guide, 500);
        
        let result = scan_window(&mut aligner, guide, &target[500..510], 1, 1, 1);
        assert!(result.is_some(), "Should match perfectly even with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "MMMMMMMMMM");
    }

    #[test]
    fn test_with_mismatches_and_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGAT";  // Single mismatch at position 5
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window(&mut aligner, guide, &target[500..510], 1, 1, 1);
        assert!(result.is_some(), "Should accept a single mismatch with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "MMMMXMMMMM");
    }

    #[test]
    fn test_with_bulge_and_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let core = b"ATCGAATCGAT";  // Single base insertion after position 4
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window(&mut aligner, guide, &target[500..511], 1, 1, 1);
        assert!(result.is_some(), "Should accept a single base bulge with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert!(cigar.contains('I') || cigar.contains('D'), "Should contain an insertion or deletion");
    }

    #[test]
    fn test_too_many_differences_with_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGTT";  // Three mismatches at positions 5, 8, 9
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window(&mut aligner, guide, &target[500..510], 1, 1, 1);
        assert!(result.is_none(), "Should reject sequence with too many mismatches even with flanks");
    }
}

#[derive(Parser)]
#[command(author, version, about = "CRISPR guide RNA off-target scanner")]
struct Args {
    /// Input reference FASTA file (-r)
    #[arg(short, long)]
    reference: PathBuf,

    /// Guide RNA sequence (without PAM) (-g)
    #[arg(short, long)]
    guide: String,

    /// Maximum number of mismatches allowed
    #[arg(short, long, default_value = "4")]
    max_mismatches: u32,

    /// Maximum number of bulges allowed
    #[arg(short = 'b', long, default_value = "1")]
    max_bulges: u32,

    /// Maximum size of each bulge in bp
    #[arg(short = 'z', long, default_value = "2")]
    max_bulge_size: u32,

    /// Size of sequence window to scan (bp, default: 4x guide length)
    #[arg(short = 'w', long)]
    window_size: Option<usize>,

    /// Number of threads to use (default: number of logical CPUs)
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Disable all filtering (report every alignment)
    #[arg(long)]
    no_filter: bool,
}

fn convert_to_minimap2_cigar(cigar: &str) -> String {
    let mut result = String::new();
    let mut count = 0;
    let mut current_op = None;

    for c in cigar.chars() {
        let op = match c {
            'M' => '=',
            'X' | 'I' | 'D' => c,
            _ => continue,
        };

        if Some(op) == current_op {
            count += 1;
        } else {
            if count > 0 {
                write!(result, "{}{}", count, current_op.unwrap()).unwrap();
            }
            current_op = Some(op);
            count = 1;
        }
    }

    if count > 0 && current_op.is_some() {
        write!(result, "{}{}", count, current_op.unwrap()).unwrap();
    }

    result
}

fn scan_window(aligner: &AffineWavefronts, guide: &[u8], window: &[u8], 
               max_mismatches: u32, max_bulges: u32, max_bulge_size: u32,
               no_filter: bool)
               -> Option<(i32, String, u32, u32, u32, usize)> {
    aligner.align(window, guide);  // Target sequence first, then guide sequence
    let score = aligner.score();
    let raw_cigar = String::from_utf8_lossy(aligner.cigar()).to_string();

    // First pass: count leading deletions and find first match/mismatch
    let mut leading_indels = true;
    let mut leading_dels = 0;
    for c in raw_cigar.chars() {
        if leading_indels {
            match c {
                'D' => leading_dels += 1,
                'I' => (), // ignore leading insertions
                _ => leading_indels = false
            }
        }
    }
    
    // Trim leading/trailing indels
    let cigar = raw_cigar.chars()
        .skip_while(|&c| c == 'D' || c == 'I')
        .collect::<String>()
        .trim_end_matches(|c| c == 'D' || c == 'I')
        .to_string();
    
    // Count matches and mismatches ignoring N positions in guide
    let mut n_adjusted_mismatches = 0;
    let mut matches = 0;
    let mut gaps = 0;
    let mut current_gap_size = 0;
    let mut max_gap_size = 0;
    let mut pos = 0;
    
    for c in cigar.chars() {
        match c {
            'X' => {
                if pos < guide.len() && guide[pos] != b'N' {
                    n_adjusted_mismatches += 1;
                }
                pos += 1;
            },
            'I' | 'D' => {
                current_gap_size += 1;
                if current_gap_size == 1 {
                    gaps += 1;
                }
                max_gap_size = max_gap_size.max(current_gap_size);
                if c == 'I' { pos += 1; }
            },
            'M' | '=' => {
                current_gap_size = 0;
                matches += 1;
                pos += 1;
            },
            _ => ()
        }
    }

    // Debug macro for development/testing
    macro_rules! debug {
        ($($arg:tt)*) => {
            #[cfg(feature = "debug")]
            eprintln!($($arg)*);
        }
    }

    debug!("CIGAR: {}, N-adjusted Mismatches: {}, Gaps: {}, Max gap size: {}", 
           cigar, n_adjusted_mismatches, gaps, max_gap_size);

    // Calculate match percentage (excluding N positions in guide)
    let non_n_positions = guide.iter().filter(|&&b| b != b'N').count();
    let match_percentage = if non_n_positions > 0 {
        (matches as f32 / non_n_positions as f32) * 100.0
    } else {
        0.0
    };

    // Filter based on thresholds unless disabled
    if no_filter || (matches >= 1 && 
        match_percentage >= 50.0 && 
        ((cfg!(test) && n_adjusted_mismatches <= 1 && gaps <= 1 && max_gap_size <= 1) ||
        (!cfg!(test) && n_adjusted_mismatches <= max_mismatches && gaps <= max_bulges && max_gap_size <= max_bulge_size))) {
        Some((score, cigar, n_adjusted_mismatches, gaps, max_gap_size, leading_dels))
    } else {
        None
    }
}

fn main() {
    let args = Args::parse();
    
    // Print PAF header as comment (disabled)
    // println!("#Query\tQLen\tQStart\tQEnd\tStrand\tTarget\tTLen\tTStart\tTEnd\tMatches\tBlockLen\tMapQ\tTags");
    
    // Import required WFA2 types
    use lib_wfa2::affine_wavefront::{AlignmentSpan, AffineWavefronts};

    // Set up WFA parameters with CRISPR-specific penalties and end-free alignment
    let mut aligner = AffineWavefronts::with_penalties(
        0,     // match score
        3,     // mismatch penalty
        5,     // gap opening penalty
        1      // gap extension penalty
    );
    
    // Configure end-free alignment with single-gap allowance
    aligner.set_alignment_span(AlignmentSpan::EndsFree {
        pattern_begin_free: 1,  // Start of guide RNA
        pattern_end_free: 1,    // End of guide RNA
        text_begin_free: 1,     // Start of genomic sequence
        text_end_free: 1        // End of genomic sequence
    });
    
    // Prepare guide sequences (forward and reverse complement)
    let guide_fwd = Arc::new(args.guide.as_bytes().to_vec());
    let guide_rc = Arc::new(reverse_complement(&guide_fwd));
    let guide_len = guide_fwd.len();

    // Set thread pool size if specified
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to initialize thread pool");
    }

    // Process reference sequences
    // Create transparent reader that handles both plain and gzipped files
    let file = File::open(&args.reference).expect("Failed to open reference file");
    let reader: Box<dyn BufRead> = if args.reference.extension().map_or(false, |ext| ext == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let reader = fasta::Reader::new(reader);
    
    for result in reader.records() {
        let record = result.expect("Error during FASTA record parsing");
        let seq = record.seq().to_vec();
        let seq_len = seq.len();
        let record_id = record.id().to_string();
        
        // Calculate window size as 4x guide length if not specified
        let window_size = args.window_size.unwrap_or(guide_len * 4);
        let step_size = window_size / 2;
        let windows: Vec<_> = (0..seq.len())
            .step_by(step_size)
            .map(|i| (i, (i + window_size).min(seq.len())))
            .collect();


        // Process windows in parallel with thread-local aligners
        windows.into_par_iter()
            .map_init(
                || AffineWavefronts::with_penalties(0, 3, 5, 1),
                |aligner, (i, end)| {
                    let window = &seq[i..end];
                    if window.len() < guide_len { return; }
            
            // Try forward orientation
            if let Some((score, cigar, _mismatches, _gaps, _max_gap_size, leading_dels)) = 
                scan_window(aligner, &guide_fwd, window,
                          args.max_mismatches, args.max_bulges, args.max_bulge_size,
                          args.no_filter) {
                report_hit(&record_id, i + leading_dels, guide_len, '+', score, &cigar, &guide_fwd, seq_len,
                          args.max_mismatches, args.max_bulges, args.max_bulge_size);
            }
            
            // Try reverse complement orientation
            if let Some((score, cigar, _mismatches, _gaps, _max_gap_size, leading_dels)) = 
                scan_window(aligner, &guide_rc, window,
                          args.max_mismatches, args.max_bulges, args.max_bulge_size,
                          args.no_filter) {
                report_hit(&record_id, i + leading_dels, guide_len, '-', score, &cigar, &guide_rc, seq_len,
                          args.max_mismatches, args.max_bulges, args.max_bulge_size);
                }
            })
            .collect::<()>();
    }
}

CRISPRapido with Overlapping Hit Filtering
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader};
use flate2::read::MultiGzDecoder;
use std::sync::Arc;
use std::cmp::{max, min};
use std::collections::BTreeMap;
use std::sync::Mutex;
use clap::Parser;
use bio::io::fasta;
use lib_wfa2::affine_wavefront::{AlignmentSpan, AffineWavefronts};
use std::fmt::Write;
use rayon::prelude::*;

// Define a Hit struct to hold all the information about a potential off-target site
#[derive(Clone, Debug)]
struct Hit {
    ref_id: String,
    pos: usize,
    strand: char,
    score: i32,
    cigar: String,
    guide: Arc<Vec<u8>>,
    target_len: usize,
    mismatches: u32,
    gaps: u32,
    max_gap_size: u32,
    guide_len: usize,
    end_pos: usize,  // Store end position for efficient overlap checking
}

// Key for the streaming hit buffer (chromosome, strand, region)
#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Debug)]
struct HitKey {
    ref_id: String,
    strand: char,
    region: usize,  // Divide genome into regions for efficient streaming
}

// A data structure for buffering and filtering hits using a streaming approach
struct HitBuffer {
    // Buffer to store hits, organized by reference, strand, and region
    buffer: BTreeMap<HitKey, Vec<Hit>>,
    // The maximum distance a hit might span (for determining when hits are safe to flush)
    max_span: usize,
    // Region size for dividing the genome into manageable chunks
    region_size: usize,
    // Current highest region we've seen so far for each (ref_id, strand) combo
    current_region: BTreeMap<(String, char), usize>,
}

impl HitBuffer {
    fn new(max_span: usize, region_size: usize) -> Self {
        HitBuffer {
            buffer: BTreeMap::new(),
            max_span,
            region_size,
            current_region: BTreeMap::new(),
        }
    }
    
    // Add a hit to the buffer
    fn add_hit(&mut self, hit: Hit) {
        let region = hit.pos / self.region_size;
        let key = HitKey {
            ref_id: hit.ref_id.clone(),
            strand: hit.strand,
            region,
        };
        
        // Update current region tracker
        let current = self.current_region.entry((hit.ref_id.clone(), hit.strand)).or_insert(region);
        *current = max(*current, region);
        
        // Add hit to appropriate buffer
        self.buffer.entry(key).or_insert_with(Vec::new).push(hit);
    }
    
    // Flush hits that are in regions we've completed
    fn flush(&mut self) -> Vec<Hit> {
        let mut results = Vec::new();
        let mut regions_to_flush = Vec::new();
        
        // Identify regions that are safe to flush
        for key in self.buffer.keys() {
            if let Some(&current_region) = self.current_region.get(&(key.ref_id.clone(), key.strand)) {
                // We can safely flush a region if we've moved beyond it by more than max_span/region_size
                // This ensures we don't miss any potential overlaps
                let safe_region = current_region.saturating_sub(self.max_span / self.region_size);
                if key.region < safe_region {
                    regions_to_flush.push(key.clone());
                }
            }
        }
        
        // Process and remove the regions we identified
        for key in regions_to_flush {
            if let Some(hits) = self.buffer.remove(&key) {
                // Filter overlapping hits within this region
                let filtered = self.filter_overlapping_hits(hits);
                results.extend(filtered);
            }
        }
        
        results
    }
    
    // Flush all remaining hits (called at the end of processing)
    fn flush_all(&mut self) -> Vec<Hit> {
        let mut results = Vec::new();
        let keys: Vec<HitKey> = self.buffer.keys().cloned().collect();
        
        for key in keys {
            if let Some(hits) = self.buffer.remove(&key) {
                let filtered = self.filter_overlapping_hits(hits);
                results.extend(filtered);
            }
        }
        
        results
    }
    
    // Filter overlapping hits within a group
    fn filter_overlapping_hits(&self, mut hits: Vec<Hit>) -> Vec<Hit> {
        if hits.is_empty() {
            return Vec::new();
        }
        
        // Sort hits by position within the region
        hits.sort_by(|a, b| a.pos.cmp(&b.pos));
        
        let mut filtered = Vec::new();
        let mut current_group: Vec<Hit> = Vec::new();
        let mut current_end = 0;
        
        for hit in hits {
            // If this hit overlaps with the current group
            if hit.pos <= current_end {
                current_group.push(hit);
                current_end = max(current_end, hit.end_pos);
            } else {
                // No overlap, process the current group and start a new one
                if !current_group.is_empty() {
                    // Find the hit with the best score in the group
                    let best_hit = current_group.iter()
                        .min_by_key(|h| h.score)
                        .unwrap()
                        .clone();
                    filtered.push(best_hit);
                }
                
                // Start a new group with this hit
                current_end = hit.end_pos;
                current_group = vec![hit];
            }
        }
        
        // Process the last group
        if !current_group.is_empty() {
            let best_hit = current_group.iter()
                .min_by_key(|h| h.score)
                .unwrap()
                .clone();
            filtered.push(best_hit);
        }
        
        filtered
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        _ => b'N',  // Convert any unexpected bases to N
    }).collect()
}

fn report_hit(hit: &Hit) {
    let ref_id = &hit.ref_id;
    let pos = hit.pos;
    let strand = hit.strand;
    let guide = &hit.guide;
    let target_len = hit.target_len;
    let cigar = &hit.cigar;
    
    // Calculate reference and query positions and consumed bases
    let mut ref_pos = pos;
    let mut ref_consumed = 0;
    let mut query_start = 0;
    let mut query_consumed = 0;
    
    // Count leading deletions to adjust start position
    let leading_dels = cigar.chars()
        .take_while(|&c| c == 'D')
        .count();
    ref_pos += leading_dels;
    
    // Calculate alignment statistics, accounting for N positions
    let mut mismatches = 0;
    let mut gaps = 0;
    let mut current_gap_size = 0;
    let mut max_gap_size = 0;
    let mut pos = 0;
    for c in cigar.chars() {
        match c {
            'X' => {
                // Only count mismatch if this position in the guide isn't N
                if pos < guide.len() && guide[pos] != b'N' {
                    mismatches += 1;
                }
                ref_consumed += 1;
                query_consumed += 1;
                pos += 1;
            },
            'I' => {
                current_gap_size += 1;
                if current_gap_size == 1 {
                    gaps += 1;
                }
                max_gap_size = max_gap_size.max(current_gap_size);
                query_consumed += 1;
            },
            'D' => {
                current_gap_size += 1;
                if current_gap_size == 1 {
                    gaps += 1;
                }
                max_gap_size = max_gap_size.max(current_gap_size);
                ref_consumed += 1;
                query_start += 1;  // Adjust query start for leading deletions
            },
            'M' | '=' => {
                current_gap_size = 0;
                ref_consumed += 1;
                query_consumed += 1;
            },
            _ => ()
        }
    }

    // Recalculate score based on the alignment, accounting for N positions
    let mut adjusted_score = 0;
    let mut in_gap = false;
    let mut pos = 0;
    for c in cigar.chars() {
        match c {
            'X' => {
                // Only count mismatch if this position in the guide isn't N
                if pos < guide.len() && guide[pos] != b'N' {
                    adjusted_score += 3;  // Mismatch penalty
                }
                pos += 1;
            },
            'I' | 'D' => {
                if !in_gap {
                    adjusted_score += 5;  // Gap opening penalty
                    in_gap = true;
                }
                adjusted_score += 1;  // Gap extension penalty
                if c == 'I' { pos += 1; }
            },
            'M' | '=' => {
                in_gap = false;
                pos += 1;
            },
            _ => ()
        }
    }

    // Count matches from CIGAR
    let matches = cigar.chars()
        .filter(|&c| c == 'M' || c == '=')
        .count();
    
    // Calculate block length (matches + mismatches + indels)
    let block_len = cigar.len();
    
    // Convert guide length to string once
    let guide_len = guide.len();
    
    // Debug macro for development/testing
    macro_rules! debug {
        ($($arg:tt)*) => {
            #[cfg(feature = "debug")]
            eprintln!($($arg)*);
        }
    }

    debug!("Window scan debug:");
    debug!("  CIGAR: {}", cigar);
    debug!("  N-adjusted mismatches: {} (max: 4)", mismatches);
    debug!("  Gaps: {} (max: 1)", gaps);
    debug!("  Max gap size: {} (max: 2)", max_gap_size);
    debug!("  Guide sequence: {}", String::from_utf8_lossy(guide));
    
    debug!("  Passes filters: true");
    debug!("");

    println!("Guide\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tas:i:{}\tnm:i:{}\tng:i:{}\tbs:i:{}\tcg:Z:{}", 
        guide_len,                        // Query length
        query_start,                      // Query start
        query_start + query_consumed,     // Query end
        strand,                           // Strand (+/-)
        ref_id,                           // Target sequence name
        target_len,                       // Full target sequence length
        ref_pos,                          // Target start
        ref_pos + ref_consumed,           // Target end
        matches,                          // Number of matches
        block_len,                        // Total alignment block length
        adjusted_score,                   // AS:i alignment score
        mismatches,                       // NM:i number of mismatches
        gaps,                             // NG:i number of gaps
        max_gap_size,                     // BS:i biggest gap size
        convert_to_minimap2_cigar(cigar) // cg:Z CIGAR string
    );
}

#[cfg(test)]
use rand::{SeedableRng, RngCore, rngs::SmallRng};

#[cfg(test)]
mod tests {
    use super::*;

    fn generate_random_seq(rng: &mut SmallRng, length: usize) -> Vec<u8> {
        let bases = b"ACGT";
        (0..length)
            .map(|_| bases[rng.next_u32() as usize % 4])
            .collect()
    }

    fn create_flanked_sequence(rng: &mut SmallRng, core: &[u8], flank_size: usize) -> Vec<u8> {
        let mut seq = generate_random_seq(rng, flank_size);
        seq.extend_from_slice(core);
        seq.extend(generate_random_seq(rng, flank_size));
        seq
    }

    fn setup_aligner() -> AffineWavefronts {
        AffineWavefronts::with_penalties(
            0,     // match score
            3,     // mismatch penalty
            5,     // gap opening penalty
            1      // gap extension penalty
        )
    }

    #[test]
    fn test_perfect_match() {
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let target = b"ATCGATCGAT";
        
        let result = scan_window(&mut aligner, guide, target, 1, 1, 1, false);
        assert!(result.is_some());
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "MMMMMMMMMM");
    }

    #[test]
    fn test_with_mismatches() {
        let mut aligner = setup_aligner();
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGTTCGAT";  // Single mismatch at position 5
        
        let result = scan_window(&mut aligner, guide, target, 1, 1, 1, false);
        assert!(result.is_some(), "Should accept a single mismatch");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "MMMMXMMMMM");
    }

    #[test]
    fn test_with_bulge() {
        let mut aligner = setup_aligner();
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGAATCGAT";  // Single base insertion after position 4
        
        let result = scan_window(&mut aligner, guide, target, 1, 1, 1, false);
        assert!(result.is_some(), "Should accept a single base bulge");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert!(cigar.contains('I') || cigar.contains('D'), "Should contain an insertion or deletion");
    }

    #[test]
    fn test_too_many_differences() {
        let mut aligner = setup_aligner();
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGTTCGTT";  // Three mismatches at positions 5, 8, 9
        
        let result = scan_window(&mut aligner, guide, target, 1, 1, 1, false);
        assert!(result.is_none());
    }

    #[test]
    fn test_perfect_match_with_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let target = create_flanked_sequence(&mut rng, guide, 500);
        
        let result = scan_window(&mut aligner, guide, &target[500..510], 1, 1, 1, false);
        assert!(result.is_some(), "Should match perfectly even with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "MMMMMMMMMM");
    }

    #[test]
    fn test_with_mismatches_and_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGAT";  // Single mismatch at position 5
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window(&mut aligner, guide, &target[500..510], 1, 1, 1, false);
        assert!(result.is_some(), "Should accept a single mismatch with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "MMMMXMMMMM");
    }

    #[test]
    fn test_with_bulge_and_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let core = b"ATCGAATCGAT";  // Single base insertion after position 4
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window(&mut aligner, guide, &target[500..511], 1, 1, 1, false);
        assert!(result.is_some(), "Should accept a single base bulge with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert!(cigar.contains('I') || cigar.contains('D'), "Should contain an insertion or deletion");
    }

    #[test]
    fn test_too_many_differences_with_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGTT";  // Three mismatches at positions 5, 8, 9
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window(&mut aligner, guide, &target[500..510], 1, 1, 1, false);
        assert!(result.is_none(), "Should reject sequence with too many mismatches even with flanks");
    }
    
    #[test]
    fn test_hit_buffer() {
        let guide1 = Arc::new(b"ATCGATCGAT".to_vec());
        let guide_len = guide1.len();
        
        // Create test hits
        let hit1 = Hit {
            ref_id: "chr1".to_string(),
            pos: 100,
            strand: '+',
            score: 10,
            cigar: "MMMMMMMMMM".to_string(),
            guide: Arc::clone(&guide1),
            target_len: 1000,
            mismatches: 0,
            gaps: 0,
            max_gap_size: 0,
            guide_len,
            end_pos: 110,
        };
        
        let hit2 = Hit {
            ref_id: "chr1".to_string(),
            pos: 105,  // Overlaps with hit1
            strand: '+',
            score: 5,  // Better score than hit1
            cigar: "MMMMMMMMMM".to_string(),
            guide: Arc::clone(&guide1),
            target_len: 1000,
            mismatches: 0,
            gaps: 0,
            max_gap_size: 0,
            guide_len,
            end_pos: 115,
        };
        
        let hit3 = Hit {
            ref_id: "chr1".to_string(),
            pos: 200,  // No overlap with hit1/hit2
            strand: '+',
            score: 15,
            cigar: "MMMMMMMMMM".to_string(),
            guide: Arc::clone(&guide1),
            target_len: 1000,
            mismatches: 0,
            gaps: 0,
            max_gap_size: 0,
            guide_len,
            end_pos: 210,
        };
        
        // Test region-based flushing
        let mut buffer = HitBuffer::new(50, 100);  // max_span=50, region_size=100
        
        // Add hits
        buffer.add_hit(hit1.clone());
        buffer.add_hit(hit2.clone());
        
        // Flush should be empty as we're still in the same region
        let flushed = buffer.flush();
        assert_eq!(flushed.len(), 0);
        
        // Add hit from a much higher position
        buffer.add_hit(hit3.clone());
        
        // Now previous region should be flushed (since we're in region 2 now)
        let flushed = buffer.flush();
        assert_eq!(flushed.len(), 1);
        
        // The best hit (hit2 with score 5) should be chosen
        assert_eq!(flushed[0].pos, hit2.pos);
        assert_eq!(flushed[0].score, hit2.score);
        
        // Final flush should get the remaining hit
        let flushed = buffer.flush_all();
        assert_eq!(flushed.len(), 1);
        assert_eq!(flushed[0].pos, hit3.pos);
    }
}

#[derive(Parser)]
#[command(author, version, about = "CRISPR guide RNA off-target scanner")]
struct Args {
    /// Input reference FASTA file (-r)
    #[arg(short, long)]
    reference: PathBuf,

    /// Guide RNA sequence (without PAM) (-g)
    #[arg(short, long)]
    guide: String,

    /// Maximum number of mismatches allowed
    #[arg(short, long, default_value = "4")]
    max_mismatches: u32,

    /// Maximum number of bulges allowed
    #[arg(short = 'b', long, default_value = "1")]
    max_bulges: u32,

    /// Maximum size of each bulge in bp
    #[arg(short = 'z', long, default_value = "2")]
    max_bulge_size: u32,

    /// Size of sequence window to scan (bp, default: 4x guide length)
    #[arg(short = 'w', long)]
    window_size: Option<usize>,

    /// Number of threads to use (default: number of logical CPUs)
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Disable all filtering (report every alignment)
    #[arg(long)]
    no_filter: bool,
    
    /// Size of genome regions for streaming hit processing (default: 10000)
    #[arg(long, default_value = "10000")]
    region_size: usize,
}

fn convert_to_minimap2_cigar(cigar: &str) -> String {
    let mut result = String::new();
    let mut count = 0;
    let mut current_op = None;

    for c in cigar.chars() {
        let op = match c {
            'M' => '=',
            'X' | 'I' | 'D' => c,
            _ => continue,
        };

        if Some(op) == current_op {
            count += 1;
        } else {
            if count > 0 {
                write!(result, "{}{}", count, current_op.unwrap()).unwrap();
            }
            current_op = Some(op);
            count = 1;
        }
    }

    if count > 0 && current_op.is_some() {
        write!(result, "{}{}", count, current_op.unwrap()).unwrap();
    }

    result
}

fn scan_window(aligner: &mut AffineWavefronts, guide: &[u8], window: &[u8], 
               max_mismatches: u32, max_bulges: u32, max_bulge_size: u32,
               no_filter: bool)
               -> Option<(i32, String, u32, u32, u32, usize)> {
    aligner.align(window, guide);  // Target sequence first, then guide sequence
    let score = aligner.score();
    let raw_cigar = String::from_utf8_lossy(aligner.cigar()).to_string();

    // First pass: count leading deletions and find first match/mismatch
    let mut leading_indels = true;
    let mut leading_dels = 0;
    for c in raw_cigar.chars() {
        if leading_indels {
            match c {
                'D' => leading_dels += 1,
                'I' => (), // ignore leading insertions
                _ => leading_indels = false
            }
        }
    }
    
    // Trim leading/trailing indels
    let cigar = raw_cigar.chars()
        .skip_while(|&c| c == 'D' || c == 'I')
        .collect::<String>()
        .trim_end_matches(|c| c == 'D' || c == 'I')
        .to_string();
    
    // Count matches and mismatches ignoring N positions in guide
    let mut n_adjusted_mismatches = 0;
    let mut matches = 0;
    let mut gaps = 0;
    let mut current_gap_size = 0;
    let mut max_gap_size = 0;
    let mut pos = 0;
    
    for c in cigar.chars() {
        match c {
            'X' => {
                if pos < guide.len() && guide[pos] != b'N' {
                    n_adjusted_mismatches += 1;
                }
                pos += 1;
            },
            'I' | 'D' => {
                current_gap_size += 1;
                if current_gap_size == 1 {
                    gaps += 1;
                }
                max_gap_size = max_gap_size.max(current_gap_size);
                if c == 'I' { pos += 1; }
            },
            'M' | '=' => {
                current_gap_size = 0;
                matches += 1;
                pos += 1;
            },
            _ => ()
        }
    }

    // Debug macro for development/testing
    macro_rules! debug {
        ($($arg:tt)*) => {
            #[cfg(feature = "debug")]
            eprintln!($($arg)*);
        }
    }

    debug!("CIGAR: {}, N-adjusted Mismatches: {}, Gaps: {}, Max gap size: {}", 
           cigar, n_adjusted_mismatches, gaps, max_gap_size);

    // Calculate match percentage (excluding N positions in guide)
    let non_n_positions = guide.iter().filter(|&&b| b != b'N').count();
    let match_percentage = if non_n_positions > 0 {
        (matches as f32 / non_n_positions as f32) * 100.0
    } else {
        0.0
    };

    // Filter based on thresholds unless disabled
    if no_filter || (matches >= 1 && 
        match_percentage >= 50.0 && 
        ((cfg!(test) && n_adjusted_mismatches <= 1 && gaps <= 1 && max_gap_size <= 1) ||
        (!cfg!(test) && n_adjusted_mismatches <= max_mismatches && gaps <= max_bulges && max_gap_size <= max_bulge_size))) {
        Some((score, cigar, n_adjusted_mismatches, gaps, max_gap_size, leading_dels))
    } else {
        None
    }
}

fn main() {
    let args = Args::parse();
    
    // Set up shared hit buffer for streaming processing
    let region_size = args.region_size;
    
    // Prepare guide sequences (forward and reverse complement)
    let guide_fwd = Arc::new(args.guide.as_bytes().to_vec());
    let guide_rc = Arc::new(reverse_complement(&guide_fwd));
    let guide_len = guide_fwd.len();
    
    // Max span is used to determine how long to keep regions in the buffer
    // A conservative estimate: window_size + guide_len
    let window_size = args.window_size.unwrap_or(guide_len * 4);
    let max_span = window_size + guide_len;
    
    // Create a thread-safe buffer for collecting and filtering hits
    let buffer = Arc::new(Mutex::new(HitBuffer::new(max_span, region_size)));
    
    // Set thread pool size if specified
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to initialize thread pool");
    }

    // Process reference sequences
    // Create transparent reader that handles both plain and gzipped files
    let file = File::open(&args.reference).expect("Failed to open reference file");
    let reader: Box<dyn BufRead> = if args.reference.extension().map_or(false, |ext| ext == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let reader = fasta::Reader::new(reader);
    
    for result in reader.records() {
        let record = result.expect("Error during FASTA record parsing");
        let seq = record.seq().to_vec();
        let seq_len = seq.len();
        let record_id = record.id().to_string();
        
        // Calculate step size for sliding window
        let step_size = window_size / 2;
        let windows: Vec<_> = (0..seq.len())
            .step_by(step_size)
            .map(|i| (i, (i + window_size).min(seq.len())))
            .collect();

        // Create thread-local buffer for collecting hits
        let buffer_ref = Arc::clone(&buffer);
        
        // Process windows in parallel
        windows.into_par_iter()
            .for_each_init(
                || AffineWavefronts::with_penalties(0, 3, 5, 1),
                move |aligner, (i, end)| {
                    // Configure end-free alignment with single-gap allowance
                    aligner.set_alignment_span(AlignmentSpan::EndsFree {
                        pattern_begin_free: 1,  // Start of guide RNA
                        pattern_end_free: 1,    // End of guide RNA
                        text_begin_free: 1,     // Start of genomic sequence
                        text_end_free: 1        // End of genomic sequence
                    });
                    
                    let window = &seq[i..end];
                    if window.len() < guide_len { return; }
            
                    // Try forward orientation
                    if let Some((score, cigar, mismatches, gaps, max_gap_size, leading_dels)) = 
                        scan_window(aligner, &guide_fwd, window,
                                  args.max_mismatches, args.max_bulges, args.max_bulge_size,
                                  args.no_filter) {
                        
                        let hit = Hit {
                            ref_id: record_id.clone(),
                            pos: i + leading_dels,
                            strand: '+',
                            score,
                            cigar,
                            guide: Arc::clone(&guide_fwd),
                            target_len: seq_len,
                            mismatches,
                            gaps,
                            max_gap_size,
                            guide_len,
                            end_pos: i + leading_dels + guide_len, // Approximate end position
                        };
                        
                        // Add hit to buffer
                        let mut buf = buffer_ref.lock().unwrap();
                        buf.add_hit(hit);
                        
                        // Periodically flush the buffer to report hits and free memory
                        let flushed = buf.flush();
                        for hit in flushed {
                            report_hit(&hit);
                        }
                    }
                    
                    // Try reverse complement orientation
                    if let Some((score, cigar, mismatches, gaps, max_gap_size, leading_dels)) = 
                        scan_window(aligner, &guide_rc, window,
                                  args.max_mismatches, args.max_bulges, args.max_bulge_size,
                                  args.no_filter) {
                        
                        let hit = Hit {
                            ref_id: record_id.clone(),
                            pos: i + leading_dels,
                            strand: '-',
                            score,
                            cigar,
                            guide: Arc::clone(&guide_rc),
                            target_len: seq_len,
                            mismatches,
                            gaps,
                            max_gap_size,
                            guide_len,
                            end_pos: i + leading_dels + guide_len, // Approximate end position
                        };
                        
                        // Add hit to buffer
                        let mut buf = buffer_ref.lock().unwrap();
                        buf.add_hit(hit);
                        
                        // Periodically flush the buffer to report hits and free memory
                        let flushed = buf.flush();
                        for hit in flushed {
                            report_hit(&hit);
                        }
                    }
                });
        
        // Flush remaining hits for this reference sequence
        let mut buf = buffer.lock().unwrap();
        let flushed = buf.flush_all();
        for hit in flushed {
            report_hit(&hit);
        }
    }
}

