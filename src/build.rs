use std::process::Command;
use std::str;

fn main() {
    let mut version = env!("CARGO_PKG_VERSION").to_string();
    match Command::new("git").args(&["describe", "--always"]).output() {
        Ok(output) => {
            let git = str::from_utf8(&output.stdout).unwrap().trim();
            version.push_str("-");
            version.push_str(git);
        },
        Err(_) => {}
    }

	println!("cargo:rustc-env=LUMOL_FULL_GIT_VERSION={}", version);
	println!("cargo:rerun-if-changed=(nonexistentfile)");
}
