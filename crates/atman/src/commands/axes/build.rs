use anyhow::Result;
use clap::Args as ClapArgs;

#[derive(ClapArgs, Debug)]
pub struct BuildArgs {}

pub fn run(_args: BuildArgs) -> Result<()> {
    anyhow::bail!("axes build is implemented in Task 4")
}
