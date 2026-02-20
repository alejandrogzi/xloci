use clap::Parser;
use log::info;
use simple_logger::init_with_level;
use xloci::{Args, xloci};

fn main() {
    let args = Args::parse();

    init_with_level(args.level).unwrap_or_else(|e| panic!("{}", e));
    info!("Starting xloci with args: {}", args);

    xloci(args);
}
