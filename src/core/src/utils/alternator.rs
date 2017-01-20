pub struct Alternator {
    step_count: usize,
    period: usize
}

impl Alternator {
    pub fn new(period: usize) -> Alternator {
        Alternator { step_count: 0, period: period }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alternator() {
        let mut alternator = Alternator::new(8);
        for i in 1..100 {
            assert_eq!(alternator.can_run(), i % 8 == 0);
        }
    }
}