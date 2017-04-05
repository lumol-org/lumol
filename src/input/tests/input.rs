// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license
extern crate test;
extern crate walkdir;
extern crate env_logger;

extern crate lumol;
extern crate lumol_input;

use std::{env, fs, io};
use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::prelude::*;

use walkdir::WalkDir;

use test::{TestDesc, TestDescAndFn, DynTestName, DynTestFn};
use test::ShouldPanic::No;

use lumol::sys::System;
use lumol_input::{InteractionsInput, Input, Error};

fn main() {
    env_logger::init().expect("could not set logger");
    let _cleanup = TestsCleanup;

    let args: Vec<_> = env::args().collect();
    let mut opts = match test::parse_opts(&args).expect("no options") {
        Ok(opts) => opts,
        Err(msg) => panic!("{:?}", msg),
    };
    opts.verbose = true;

    let tests = all_tests();
    let result = test::run_tests_console(&opts, tests);
    match result {
        Ok(true) => {}
        Ok(false) => std::process::exit(-1),
        Err(err) => panic!("io error when running tests: {:?}", err),
    }
}

fn all_tests() -> Vec<TestDescAndFn> {
    let mut tests = Vec::new();

    tests.extend(generate_tests("simulation/good", |path| {
        Box::new(move || {
            let input = Input::new(path.clone()).unwrap();
            input.read().unwrap();
        })
    }).expect("Could not generate the tests"));

    tests.extend(generate_tests("simulation/bad", |path| {
        Box::new(move || {
            let message = get_error_message(&path);
            let result = Input::new(path.clone())
                               .and_then(|input| input.read());

            match result {
                Err(Error::Config(reason)) => assert_eq!(reason, message),
                _ => panic!("This test should fail with a Config error")
            }
        })
    }).expect("Could not generate the tests"));

    tests.extend(generate_tests("interactions/good", |path| {
        Box::new(move || {
            let mut system = System::new();
            let input = InteractionsInput::new(path.clone()).unwrap();
            input.read(&mut system).unwrap();
        })
    }).expect("Could not generate the tests"));

    tests.extend(generate_tests("interactions/bad", |path| {
        Box::new(move || {
            let message = get_error_message(&path);

            let mut system = System::new();
            let result = InteractionsInput::new(path.clone())
                                           .and_then(|input| input.read(&mut system));

            match result {
                Err(Error::Config(reason)) => assert_eq!(reason, message),
                _ => panic!("This test should fail with a Config error")
            }
        })
    }).expect("Could not generate the tests"));

    return tests;
}

/// Generate the tests by calling `callback` for every TOML files at the given
/// `root`.
fn generate_tests<F>(root: &str, callback: F) -> Result<Vec<TestDescAndFn>, io::Error>
    where F: Fn(PathBuf) -> Box<FnMut() + Send> {
    let mut tests = Vec::new();

    let dir = PathBuf::new().join(env!("CARGO_MANIFEST_DIR"))
                            .join("tests")
                            .join(root);
    for entry in WalkDir::new(dir) {
        let entry = try!(entry);
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

                    let test = TestDescAndFn {
                        desc: TestDesc {
                            name: DynTestName(name),
                            ignore: false,
                            should_panic: No,
                        },
                        testfn: DynTestFn(callback(path.to_path_buf())),
                    };

                    tests.push(test);
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
        const REMOVE: &'static [&'static str] = &[
            "energy.dat", "filename.xyz", "cell.dat", "properties.dat",
            "file.log"
        ];

        for file in REMOVE {
            if let Err(err) = fs::remove_file(file) {
                match err.kind() {
                    io::ErrorKind::NotFound => {},
                    _ => panic!("io error in cleanup code: {}", err)
                }
            }
        }
    }
}

fn get_error_message(path: &Path) -> String {
    let mut buffer = String::new();
    File::open(path)
          .and_then(|mut file| file.read_to_string(&mut buffer))
          .expect("Could not read the input file");

    for line in buffer.lines() {
        let line = line.trim();
        if line.starts_with("#^ ") {
            return String::from(&line[3..])
        }
    }

    panic!("No error message found in {}. Please add one with the '#^ <message>' syntax.", path.display());
}
