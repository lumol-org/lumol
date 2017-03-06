// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

extern crate test;
extern crate walkdir;

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

    let args: Vec<_> = env::args().collect();
    let mut opts = match test::parse_opts(&args) {
        Some(Ok(opts)) => opts,
        Some(Err(msg)) => panic!("{:?}", msg),
        None => return,
    };
    opts.verbose = true;

    let result = test::run_tests_console(&opts, tests);
    cleanup();
    match result {
        Ok(true) => {}
        Ok(false) => std::process::exit(-1),
        Err(err) => panic!("io error when running tests: {:?}", err),
    }
}

fn generate_tests<F>(group: &str, callback: F) -> Result<Vec<TestDescAndFn>, io::Error>
    where F: Fn(PathBuf) -> Box<FnMut() + Send> {
    let mut tests = Vec::new();

    let dir = PathBuf::new().join(env!("CARGO_MANIFEST_DIR")).join("tests").join(group);
    for entry in WalkDir::new(dir) {
        let entry = try!(entry);
        let file_type = entry.file_type();
        if file_type.is_file() {
            if let Some(extension) = entry.path().extension() {
                if extension == "toml" {
                    let path = entry.path();
                    let name = String::from(group) + "/" + path.file_name().expect("Missing file name")
                                                               .to_str().expect("File name is invalid UTF-8");

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

fn cleanup() {
    const REMOVE: &'static [&'static str] = &[
        "energy.dat", "filename.xyz", "cell.dat", "properties.dat"
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
