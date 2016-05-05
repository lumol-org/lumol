// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

/// Extract the table at the given `key`, from a TOML table interpreted as
/// a `context`
///
/// ```
/// let propagator = extract_table!("propagator", config as "simulation");
/// let control = extract_table!("control", config as "molecular dynamics");
/// ```
macro_rules! extract_table {
    ($key: expr, $config: ident as $context: expr) => ({
        let table = try!($config.get($key).ok_or(Error::from(
            format!("Missing '{}' key in {}", $key, $context)
        )));
        try!(table.as_table().ok_or(Error::from(
            format!("'{}' must be a table in {}", $key, $context)
        )))
    });
}

/// Extract the string at the given `key`, from a TOML table interpreted as
/// a `context`
macro_rules! extract_str {
    ($key: expr, $config: ident as $context: expr) => ({
        let string = try!($config.get($key).ok_or(Error::from(
            format!("Missing '{}' key in {}", $key, $context)
        )));
        try!(string.as_str().ok_or(Error::from(
            format!("'{}' must be a string in {}", $key, $context)
        )))
    });
}

/// Extract a number (integer or float) at the given `key`, from a TOML table
/// interpreted as a `context`
macro_rules! extract_number {
    ($key: expr, $config: ident as $context: expr) => ({
        let number = try!($config.get($key).ok_or(Error::from(
            format!("Missing '{}' key in {}", $key, $context)
        )));
        match *number {
            ::toml::Value::Integer(v) => v as f64,
            ::toml::Value::Float(v) => v,
            _ => return Err(Error::from(
                format!("'{}' must be a string in {}", $key, $context)
            ))
        }
    });
}

/// Extract the string 'type' key in a TOML table
///
/// ```
/// let value = match extract_type!(config) {
///     "foo" => create_foo(config),
///     "bar" => create_bar(config),
///     other => error!("Unknown type of stuff '{}'", other),
/// }
/// ```
macro_rules! extract_type {
    ($config: ident) => ({
        let typ = try!($config.get("type").ok_or(Error::from(
            format!("Missing 'type' key in {}", stringify!($config))
        )));
        try!(typ.as_str().ok_or(Error::from(
            format!("'type' key must be a string in {}", stringify!($config))
        )))
    });
}
