// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use toml::value::{Table, Value};

use crate::error::Error;

/// Extract the table at the given `key`, from the `config` TOML table
/// interpreted as a `context`.
pub fn table<'a>(key: &str, config: &'a Table, context: &str) -> Result<&'a Table, Error> {
    let table = config.get(key).ok_or(
        Error::from(format!("missing '{key}' key in {context}"))
    )?;
    return table.as_table().ok_or(
        Error::from(format!("'{key}' must be a table in {context}"))
    );
}

/// Extract the string at the given `key`, from the `config` TOML table
/// interpreted as a `context`
pub fn str<'a>(key: &str, config: &'a Table, context: &str) -> Result<&'a str, Error> {
    let string = config.get(key).ok_or(
        Error::from(format!("missing '{key}' key in {context}"))
    )?;
    return string.as_str().ok_or(
        Error::from(format!("'{key}' must be a string in {context}"))
    );
}

/// Extract a number (integer or float) at the given `key`, from the `config`
/// TOML table interpreted as a `context`
pub fn number(key: &str, config: &Table, context: &str) -> Result<f64, Error> {
    let number = config.get(key).ok_or(
        Error::from(format!("missing '{key}' key in {context}"))
    )?;
    match *number {
        ::toml::Value::Integer(v) => Ok(v as f64),
        ::toml::Value::Float(v) => Ok(v),
        _ => Err(Error::from(format!("'{key}' must be a number in {context}"))),
    }
}

/// Extract a unsigned integer at the given `key`, from the `config`
/// TOML table interpreted as a `context`
pub fn uint(key: &str, config: &Table, context: &str) -> Result<u64, Error> {
    let number = config.get(key).ok_or(
        Error::from(format!("missing '{key}' key in {context}"))
    )?;
    match *number {
        ::toml::Value::Integer(v) => {
            if v < 0 {
                Err(Error::from(format!("'{key}' must be a positive integer in {context}")))
            } else {
                Ok(v as u64)
            }
        }
        _ => Err(Error::from(format!("'{key}' must be a positive integer in {context}"))),
    }
}

/// Extract an array at the given `key`, from the `config` TOML table
/// interpreted as a `context`
pub fn slice<'a>(key: &str, config: &'a Table, context: &str) -> Result<&'a [Value], Error> {
    let array = config.get(key).ok_or(
        Error::from(format!("missing '{key}' key in {context}"))
    )?;
    let array = array.as_array().ok_or(
        Error::from(format!("'{key}' must be an array in {context}"))
    );
    return array.map(|arr| arr.as_slice());
}

/// Extract the string 'type' key in a TOML table
pub fn typ<'a>(config: &'a Table, context: &str) -> Result<&'a str, Error> {
    let typ = config.get("type").ok_or(
        Error::from(format!("missing 'type' key in {context}"))
    )?;
    return typ.as_str().ok_or(Error::from(format!("'type' key must be a string in {context}")));
}
