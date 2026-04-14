use anyhow::Result;
use clap::Args as ClapArgs;

#[derive(ClapArgs, Debug)]
pub struct Args {}

pub fn run(_args: Args) -> Result<()> {
    anyhow::bail!("qc: not yet implemented (Task 15)")
}
