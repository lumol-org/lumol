extern crate test;
extern crate walkdir;

extern crate lumol;
extern crate lumol_input;

use std::{env, fs, io};
use std::path::PathBuf;
use walkdir::WalkDir;

use test::{TestDesc, TestDescAndFn, DynTestName, DynTestFn};
use test::ShouldPanic::No;

use lumol::sys::System;
use lumol_input::{InteractionsInput, Input};

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
            assert!(Input::new(path.clone()).and_then(|input| input.read()).is_err());
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
            let mut system = System::new();
            assert!(InteractionsInput::new(path.clone()).and_then(|input| input.read(&mut system)).is_err());
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

    return Ok(tests);
}

fn cleanup() {
    const REMOVE: &'static [&'static str] = &["energy.dat", "filename.xyz"];

    for file in REMOVE {
        if let Err(err) = fs::remove_file(file) {
            match err.kind() {
                io::ErrorKind::NotFound => {},
                _ => panic!("io error in cleanup code: {}", err)
            }
        }
    }
}
