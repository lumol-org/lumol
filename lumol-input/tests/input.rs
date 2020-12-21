// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
#![allow(clippy::needless_return)]

use std::{env, fs, io};
use std::fs::File;
use std::io::prelude::*;
use std::path::PathBuf;

use walkdir::WalkDir;

use rustc_test::{DynTestFn, DynTestName, TestDesc, TestDescAndFn};
use rustc_test::ShouldPanic::No;

use lumol_core::System;
use lumol_input::{Error, Input, InteractionsInput};

fn main() {
    env_logger::init();
    let _cleanup = TestsCleanup;

    let args: Vec<_> = env::args()
                           .filter(|arg| !arg.contains("test-threads"))
                           .collect();

    let mut opts = match rustc_test::parse_opts(&args).expect("no options") {
        Ok(opts) => opts,
        Err(msg) => panic!("{:?}", msg),
    };
    opts.verbose = true;

    let tests = all_tests();
    let result = rustc_test::run_tests_console(&opts, tests);
    match result {
        Ok(true) => {}
        Ok(false) => std::process::exit(-1),
        Err(err) => panic!("io error when running tests: {:?}", err),
    }
}

fn all_tests() -> Vec<TestDescAndFn> {
    let mut tests = Vec::new();

    tests.extend(
        generate_tests("simulation/good", |path, content| {
            Box::new(move || {
                let input = Input::from_str(path.clone(), &content).unwrap();
                input.read().unwrap();
            })
        }).expect("Could not generate the tests"),
    );

    tests.extend(
        generate_tests("simulation/bad", |path, content| {
            Box::new(move || {
                let message = get_error_message(&content);
                let result = Input::from_str(path.clone(), &content).and_then(|input| input.read());

                match result {
                    Err(Error::Config(reason)) => assert_eq!(reason, message),
                    _ => panic!("This test should fail with a Config error"),
                }
            })
        }).expect("Could not generate the tests"),
    );

    tests.extend(
        generate_tests("interactions/good", |_, content| {
            Box::new(move || {
                let mut system = System::new();
                let input = InteractionsInput::from_str(&content).unwrap();
                input.read(&mut system).unwrap();
            })
        }).expect("Could not generate the tests"),
    );

    tests.extend(
        generate_tests("interactions/bad", |_, content| {
            Box::new(move || {
                let message = get_error_message(&content);

                let mut system = System::new();
                let result = InteractionsInput::from_str(&content).and_then(|input| input.read(&mut system));

                match result {
                    Err(Error::Config(reason)) => assert_eq!(reason, message),
                    _ => panic!("This test should fail with a Config error"),
                }
            })
        }).expect("Could not generate the tests"),
    );

    return tests;
}

/// Generate the tests by calling `callback` for every TOML files at the given
/// `root`.
fn generate_tests<F>(root: &str, callback: F) -> Result<Vec<TestDescAndFn>, io::Error>
where
    F: Fn(PathBuf, String) -> Box<dyn FnMut() + Send>,
{
    let mut tests = Vec::new();

    let dir = PathBuf::new().join(env!("CARGO_MANIFEST_DIR")).join("tests").join(root);
    for entry in WalkDir::new(dir) {
        let entry = entry?;
        let file_type = entry.file_type();
        if file_type.is_file() {
            if let Some(extension) = entry.path().extension() {
                if extension == "toml" {
                    let path = entry.path();
                    let name = String::from(root) + "/";
                    let name = name + path.file_name()
                                          .expect("Missing file name")
                                          .to_str()
                                          .expect("File name is invalid UTF-8");

                    let mut content = String::new();
                    File::open(path).and_then(|mut file| file.read_to_string(&mut content))
                                    .expect("Could not read the input file");

                    let count = content.split("+++").count();
                    for (i, test_case) in content.split("+++").enumerate() {
                        let name = if count > 1 {
                            format!("{} - {}/{}", name, i + 1, count)
                        } else {
                            name.clone()
                        };
                        let test = TestDescAndFn {
                            desc: TestDesc {
                                name: DynTestName(name),
                                ignore: false,
                                should_panic: No,
                                allow_fail: false
                            },
                            testfn: DynTestFn(callback(path.to_path_buf(), test_case.into())),
                        };
                        tests.push(test);
                    }
                }
            }
        }
    }

    Ok(tests)
}

/// Cleanup temporary files after the tests
struct TestsCleanup;
impl Drop for TestsCleanup {
    fn drop(&mut self) {
        const REMOVE: &[&str] = &[
            "energy.dat",
            "filename.xyz",
            "cell.dat",
            "properties.dat",
            "file.log",
            "custom.dat",
            "stress.dat",
            "forces.xyz",
        ];

        for file in REMOVE {
            if let Err(err) = fs::remove_file(file) {
                match err.kind() {
                    io::ErrorKind::NotFound => {}
                    _ => panic!("io error in cleanup code: {}", err),
                }
            }
        }
    }
}

fn get_error_message(content: &str) -> String {
    for line in content.lines() {
        let line = line.trim();
        if let Some(message) = line.strip_prefix("#^ ") {
            return message.into();
        }
    }

    panic!("No error message found. Please add one with the '#^ <message>' syntax.");
}
