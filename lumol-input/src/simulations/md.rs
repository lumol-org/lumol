// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
#![allow(clippy::wildcard_imports)]

use toml::value::Table;

use lumol_sim::md::*;
use lumol_core::units;

use crate::alternator::Alternator;
use crate::{Error, FromToml, FromTomlWithData};
use crate::extract;

impl FromToml for MolecularDynamics {
    fn from_toml(config: &Table) -> Result<MolecularDynamics, Error> {
        // Get the timestep of the simulation
        let timestep = extract::str("timestep", config, "molecular dynamics propagator")?;
        let timestep = units::from_str(timestep)?;

        let mut md;
        if let Some(integrator) = config.get("integrator") {
            let integrator = integrator.as_table().ok_or(
                Error::from("'integrator' must be a table in molecular dynamics")
            )?;

            let integrator: Box<dyn Integrator> = match extract::typ(integrator, "integrator")? {
                "BerendsenBarostat" => {
                    Box::new(BerendsenBarostat::from_toml(integrator, timestep)?)
                }
                "AnisoBerendsenBarostat" => {
                    Box::new(AnisoBerendsenBarostat::from_toml(integrator, timestep)?)
                }
                "Verlet" => Box::new(Verlet::from_toml(integrator, timestep)?),
                "VelocityVerlet" => Box::new(VelocityVerlet::from_toml(integrator, timestep)?),
                "LeapFrog" => Box::new(LeapFrog::from_toml(integrator, timestep)?),
                other => return Err(Error::from(format!("unknown integrator '{other}'"))),
            };

            md = MolecularDynamics::from_integrator(integrator);
        } else {
            md = MolecularDynamics::new(timestep);
        }

        if let Some(thermostat) = config.get("thermostat") {
            let thermostat = thermostat.as_table().ok_or(
                Error::from("'thermostat' must be a table in molecular dynamics")
            )?;

            let thermostat: Box<dyn Thermostat> = match extract::typ(thermostat, "thermostat")? {
                "Berendsen" => Box::new(BerendsenThermostat::from_toml(thermostat)?),
                "Rescale" => Box::new(RescaleThermostat::from_toml(thermostat)?),
                "CSVR" => Box::new(CSVRThermostat::from_toml(thermostat)?),
                other => return Err(Error::from(format!("unknown thermostat type '{other}'"))),
            };
            md.set_thermostat(thermostat);
        }

        if let Some(controls) = config.get("controls") {
            let controls = controls.as_array().ok_or(
                Error::from("'controls' must be an array of tables in molecular dynamics")
            )?;

            for control in controls {
                let control = control.as_table().ok_or(
                    Error::from("'controls' must be an array of tables in molecular dynamics")
                )?;

                let control: Box<dyn Control> = match extract::typ(control, "control")? {
                    "RemoveTranslation" => {
                        Box::new(Alternator::<RemoveTranslation>::from_toml(control)?)
                    }
                    "RemoveRotation" => {
                        Box::new(Alternator::<RemoveRotation>::from_toml(control)?)
                    }
                    "Rewrap" => Box::new(Alternator::<Rewrap>::from_toml(control)?),
                    other => return Err(Error::from(format!("unknown control '{other}'"))),
                };
                md.add_control(control);
            }
        }

        Ok(md)
    }
}

impl FromTomlWithData for Verlet {
    type Data = f64;
    fn from_toml(_: &Table, timestep: f64) -> Result<Verlet, Error> {
        Ok(Verlet::new(timestep))
    }
}

impl FromTomlWithData for VelocityVerlet {
    type Data = f64;
    fn from_toml(_: &Table, timestep: f64) -> Result<VelocityVerlet, Error> {
        Ok(VelocityVerlet::new(timestep))
    }
}

impl FromTomlWithData for LeapFrog {
    type Data = f64;
    fn from_toml(_: &Table, timestep: f64) -> Result<LeapFrog, Error> {
        Ok(LeapFrog::new(timestep))
    }
}

impl FromTomlWithData for BerendsenBarostat {
    type Data = f64;
    fn from_toml(config: &Table, timestep: f64) -> Result<BerendsenBarostat, Error> {
        let pressure = extract::str("pressure", config, "Berendsen barostat")?;
        let pressure = units::from_str(pressure)?;
        let tau = extract::number("timestep", config, "Berendsen barostat")?;
        Ok(BerendsenBarostat::new(timestep, pressure, tau))
    }
}

impl FromTomlWithData for AnisoBerendsenBarostat {
    type Data = f64;
    fn from_toml(config: &Table, timestep: f64) -> Result<AnisoBerendsenBarostat, Error> {
        let pressure = extract::str("pressure", config, "anisotropic Berendsen barostat")?;
        let pressure = units::from_str(pressure)?;
        let tau = extract::number("timestep", config, "anisotropic Berendsen barostat")?;
        Ok(AnisoBerendsenBarostat::hydrostatic(timestep, pressure, tau))
    }
}

impl FromToml for BerendsenThermostat {
    fn from_toml(config: &Table) -> Result<BerendsenThermostat, Error> {
        let temperature = extract::str("temperature", config, "Berendsen thermostat")?;
        let temperature = units::from_str(temperature)?;
        let tau = extract::number("timestep", config, "Berendsen thermostat")?;
        Ok(BerendsenThermostat::new(temperature, tau))
    }
}

impl FromToml for RescaleThermostat {
    fn from_toml(config: &Table) -> Result<RescaleThermostat, Error> {
        let temperature = extract::str("temperature", config, "rescale thermostat")?;
        let temperature = units::from_str(temperature)?;

        if let Some(tolerance) = config.get("tolerance") {
            let tolerance = tolerance.as_str().ok_or(
                Error::from("'tolerance' must be a string rescale thermostat")
            )?;
            let tolerance = units::from_str(tolerance)?;

            Ok(RescaleThermostat::with_tolerance(temperature, tolerance))
        } else {
            Ok(RescaleThermostat::new(temperature))
        }
    }
}

impl FromToml for CSVRThermostat {
    fn from_toml(config: &Table) -> Result<CSVRThermostat, Error> {
        let temperature = extract::str("temperature", config, "CSVR thermostat")?;
        let temperature = units::from_str(temperature)?;
        let tau = extract::number("timestep", config, "CSVR thermostat")?;
        Ok(CSVRThermostat::new(temperature, tau))
    }
}

impl FromToml for Alternator<RemoveTranslation> {
    fn from_toml(config: &Table) -> Result<Alternator<RemoveTranslation>, Error> {
        let every = if config.contains_key("every") {
            extract::uint("every", config, "RemoveTranslation control")?
        } else {
            1
        };
        Ok(Alternator::new(every, RemoveTranslation))
    }
}

impl FromToml for Alternator<RemoveRotation> {
    fn from_toml(config: &Table) -> Result<Alternator<RemoveRotation>, Error> {
        let every = if config.contains_key("every") {
            extract::uint("every", config, "RemoveRotation control")?
        } else {
            1
        };
        Ok(Alternator::new(every, RemoveRotation))
    }
}

impl FromToml for Alternator<Rewrap> {
    fn from_toml(config: &Table) -> Result<Alternator<Rewrap>, Error> {
        let every = if config.contains_key("every") {
            extract::uint("every", config, "Rewrap control")?
        } else {
            1
        };
        Ok(Alternator::new(every, Rewrap))
    }
}
