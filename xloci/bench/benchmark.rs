use clap::Parser;
use std::process::{Command, ExitStatus, Stdio};

const STDOUT: &str = "output.fa";
const STDOUT2: &str = "transcripts.bed";

const FA_BED_TOOLS: [&str; 4] = [
    "xloci -s {assembly} -o output -r {regions}",
    "bedtools getfasta -fi {assembly} -bed {regions} -split -name -fo output.fa",
    "gffread -w output.fa -g {assembly} --in-bed {regions}",
    "bed2gtf -b {regions} -o transcripts.gtf -n && agat_sp_extract_sequences.pl -g transcripts.gtf -f {assembly} -t exon --merge --cpu 16",
];

const FA_GTF_TOOLS: [&str; 5] = [
    "xloci -s {assembly} -o output -r {regions}",
    "gffread {regions} --bed -o transcripts.bed && cat -p --color never transcripts.bed | choose :11 -o '\\t' > tmp.bed && bedtools getfasta -fi {assembly} -bed tmp.bed -split -name -fo output.fa",
    "gxf2bed -i {regions} -o transcripts.bed && bedtools getfasta -fi {assembly} -bed transcripts.bed -split -name -fo output.fa",
    "gffread -w output.fa -g {assembly} {regions}",
    // "agat_sp_extract_sequences.pl -g {regions} -f {assembly} -t exon --merge",
    "gtf_to_fasta {regions} {assembly} > output.fa",
];

const TWOBIT_BED_TOOLS: [&str; 2] = [
    "xloci -s {assembly} -o output -r {regions}",
    "twoBitToFa -bed={regions} {assembly} output.fa",
];

const TWOBIT_GTF_TOOLS: [&str; 5] = [
    "xloci -s {assembly} -o output -r {regions}",
    "gxf2bed -i {regions} -o transcripts.bed && xloci -s {assembly} -o output -r transcripts.bed",
    "gffread {regions} --bed -o transcripts.bed && cat -p --color never transcripts.bed | choose :11  -o '\\t' > tmp.bed && xloci -s {assembly} -o output -r tmp.bed",
    "gffread {regions} --bed -o transcripts.bed && cat -p --color never transcripts.bed | choose :11  -o '\\t' > tmp.bed && twoBitToFa -bed=tmp.bed {assembly} output.fa",
    "gxf2bed -i {regions} -o transcripts.bed && sort -k1,1 -k2,2n -k3,3n transcripts.bed > tmp.bed && twoBitToFa -bed=tmp.bed {assembly} output.fa",
];

/// Command line arguments for the benchmark utility.
///
/// This struct defines configuration options for running performance benchmarks
/// of various chromosome size extraction tools using hyperfine.
#[derive(Debug, Parser)]
pub struct Args {
    /// Path to the reference directory containing test assemblies
    #[clap(short = 'f', long = "fasta", help = "Path to the assembly file")]
    fasta: String,

    #[clap(short = 't', long = "twobit", help = "Path to 2bit file")]
    twobit: String,

    #[clap(short = 'b', long = "bed", help = "Path to bed file")]
    bed: String,

    #[clap(short = 'g', long = "gtf", help = "Path to gtf file")]
    gtf: String,

    /// Additional arguments to pass to the hyperfine benchmarking tool
    #[clap(short = 'a',
        value_delimiter = ',',
        num_args = 1..,
        help = "Extra arguments to pass to hyperfine"
    )]
    hyperfine_args: Vec<String>,
}

/// Configuration for hyperfine benchmark execution.
///
/// This struct encapsulates all parameters needed to run hyperfine benchmarks,
/// including warmup runs, execution limits, output formats, and commands to test.
pub struct HyperfineCall {
    /// Number of warmup runs before actual benchmarking
    pub warmup: u32,
    /// Minimum number of benchmark runs
    pub min_runs: u32,
    /// Maximum number of benchmark runs (optional)
    pub max_runs: Option<u32>,
    /// Path to export results in CSV format
    pub export_csv: Option<String>,
    /// Path to export results in Markdown format
    pub export_markdown: Option<String>,
    /// Parameterized variables for command substitution
    pub parameters: Vec<(String, Vec<String>)>,
    /// Setup command to run before each benchmark
    pub setup: Option<String>,
    /// Cleanup command to run after each benchmark
    pub cleanup: Option<String>,
    /// List of commands to benchmark
    pub commands: Vec<String>,
    /// Additional hyperfine command line arguments
    pub extras: Vec<String>,
}

impl Default for HyperfineCall {
    /// Creates a default HyperfineCall with sensible baseline settings.
    ///
    /// Sets up reasonable defaults for benchmarking with 3 warmup runs
    /// and 5 minimum runs, suitable for most performance testing scenarios.
    fn default() -> Self {
        Self {
            warmup: 3,
            min_runs: 5,
            max_runs: None,
            export_csv: None,
            export_markdown: None,
            parameters: Vec::new(),
            setup: None,
            cleanup: None,
            commands: Vec::new(),
            extras: Vec::new(),
        }
    }
}

impl HyperfineCall {
    /// Executes the hyperfine benchmark with the configured parameters.
    ///
    /// This method builds and executes a hyperfine command using the struct's
    /// configuration. It sets up all command line arguments including warmup,
    /// runs, exports, parameters, setup/cleanup commands, and the actual
    /// benchmark commands.
    ///
    /// # Returns
    ///
    /// ExitStatus from the hyperfine process execution
    ///
    /// # Panics
    ///
    /// Panics if hyperfine command cannot be executed
    pub fn invoke(&self) -> ExitStatus {
        let mut command = Command::new("hyperfine");

        command
            .stdout(Stdio::inherit())
            .stderr(Stdio::inherit())
            .stdin(Stdio::null());

        command.arg("--warmup").arg(self.warmup.to_string());
        command.arg("--min-runs").arg(self.min_runs.to_string());
        if let Some(export_csv) = &self.export_csv {
            command.arg("--export-csv").arg(export_csv);
        }
        if let Some(export_markdown) = &self.export_markdown {
            command.arg("--export-markdown").arg(export_markdown);
        }
        for (flag, values) in &self.parameters {
            command.arg("-L").arg(flag).arg(values.join(","));
        }
        if let Some(setup) = &self.setup {
            command.arg("--setup").arg(setup);
        }
        if let Some(cleanup) = &self.cleanup {
            command.arg("--cleanup").arg(cleanup);
        }
        if let Some(max_runs) = self.max_runs {
            command.arg("--max-runs").arg(max_runs.to_string());
        }
        if !self.extras.is_empty() {
            command.args(&self.extras);
        }

        for cmd in &self.commands {
            command.arg(cmd);
        }

        command.status().expect("Failed to run hyperfine")
    }
}

/// Runs comprehensive benchmarks for chromosome size extraction tools.
///
/// This function sets up and executes hyperfine benchmarks comparing multiple
/// tools across various genome assemblies. It configures warmup runs, output formats,
/// parameterized assembly testing, and cleanup operations.
///
/// # Returns
///
/// Tuple containing (csv_path, markdown_path) for the benchmark results
///
/// # Errors
///
/// Returns error if benchmark fails or if file system operations fail
///
/// # Examples
///
/// ```ignore
/// let (csv, md) = benchmark()?;
/// println!("Results: {} {}", csv, md);
/// ```
fn benchmark() {
    let args = Args::parse();

    std::fs::create_dir_all("runs").unwrap_or_else(|e| panic!("{}", e));
    let triplets = vec![
        (&args.fasta, &args.bed, FA_BED_TOOLS.to_vec(), "fa_bed"),
        (&args.fasta, &args.gtf, FA_GTF_TOOLS.to_vec(), "fa_gtf"),
        (
            &args.twobit,
            &args.bed,
            TWOBIT_BED_TOOLS.to_vec(),
            "2bit_bed",
        ),
        (
            &args.twobit,
            &args.gtf,
            TWOBIT_GTF_TOOLS.to_vec(),
            "2bit_gtf",
        ),
    ];

    for (assembly, regions, tools, run_name) in triplets {
        let csv = format!("bench_{}.csv", run_name);
        let md = format!("bench_{}.md", run_name);

        #[allow(clippy::needless_update)]
        let code = HyperfineCall {
            warmup: 3,
            min_runs: 3,
            max_runs: Some(10),
            export_csv: Some(format!("runs/{}", csv).to_string()),
            export_markdown: Some(format!("runs/{}", md).to_string()),
            parameters: vec![
                ("assembly".to_string(), vec![assembly.to_string()]),
                ("regions".to_string(), vec![regions.to_string()]),
            ],
            setup: Some("cargo build --release".to_string()),
            cleanup: Some(format!("rm -rf output tmp.bed {} {}", STDOUT, STDOUT2)),
            commands: tools
                .iter()
                .map(|cmd| cmd.to_string())
                .collect::<Vec<String>>(),
            extras: args
                .hyperfine_args
                .iter()
                .map(|s| format!("--{}", s))
                .collect(),
            ..Default::default()
        }
        .invoke()
        .code()
        .expect("Benchmark terminated unexpectedly");

        if code != 0 {
            eprintln!("Benchmark failed with exit code {}", code);
        }
    }
}

/// Main entry point for the benchmark utility.
///
/// This function parses command line arguments, runs the benchmark suite,
/// and reports the location of result files or any errors that occurred.
///
/// # Examples
///
/// ```ignore
/// // Run with custom asset directory
/// ./chromsize-benchmark -d /path/to/assemblies
/// ```
fn main() {
    benchmark();
    // match benchmark() {
    //     Ok((csv, md)) => {
    //         println!("Benchmark results saved to:");
    //         println!("  - {}", csv);
    //         println!("  - {}", md);
    //     }
    //     Err(e) => eprintln!("{}", e),
    // }
}
