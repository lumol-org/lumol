pub struct Alternator {
    step_count: usize,
    period: usize
}

impl Alternator {
    pub fn from_frequency(frequency: f64) -> Alternator {
        if frequency <= 0.0 && frequency > 1.0 {
            panic!("Frequency must be between 0 and 1.")
        }
        Alternator { step_count: 0, period: (1.0 / frequency) as usize }
    }

    pub fn can_run(&mut self) -> bool {
        self.step_count += 1;
        if self.step_count < self.period {
            return false;
        }
        self.step_count = 0;
        true
    }
}