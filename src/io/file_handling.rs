use std::env;
use std::fs;
use is_prime::is_prime;

pub struct Config {
    file_path: String,
}

impl Config {
    pub fn new() -> Config {
        let args: Vec<String> = env::args().collect();
        let file_path = args[1].clone();
        Config { file_path }
    }
}

/// Read in msolve input files of the following format
/// 1. The **first line** consists of the names of the variables, commata separated.
/// 2. The **second line** consists of a nonnegative 32-bit integer representing the characteristic 
///    of the field.
/// 3. The **remaining lines** consist of the polynomials generating the system to be solved. The 
///    polynomials must be commata separated, no shorthand notation, i.e. use `*` between 
///    variables, coefficients etc. and use `^` for exponents of variables.
pub fn read_file(config: Config) -> (Vec<String>, u32, Vec<usize>, Vec<i32>, Vec<i32>) {
    let contents  = fs::read_to_string(config.file_path)
        .expect("Should have been able to read the file")
        .replace(" ", "");
    let lines: Vec<&str> = contents.lines().collect();
    let variables: Vec<String> = lines[0].split(",").map(String::from).collect();
    let characteristic: u32 = lines[1].parse::<u32>().unwrap_or(u32::MAX);
    if characteristic != 0 && !is_prime(&characteristic.to_string()) {
        panic!("Wrong characteristic!");
    }
    // get correct polynomials, i.e. several polys per line, one poly over several lines, etc.
    let tmp = &lines[2..lines.len()].join("");
    let polynomials: Vec<&str> = tmp.split(",").collect();
    // read in each equation
    let mut exponents: Vec<i32>    = Vec::new();
    let mut coefficients: Vec<i32> = Vec::new();
    let mut lengths: Vec<usize>    = Vec::new();
    for i in 2..lines.len() {
        lengths.push(lines[i].chars().filter(|v| *v == '+' || *v == '-').count() + 1);
        for positive_start in lines[i].split('+') {
            let monomials = positive_start.split('-').filter(|v| !v.is_empty());
            let mut sign: i32 = 1;
            for mons in monomials {
                // first entry is the only positive one
                let mon: Vec<&str> = mons.split('*').filter(|v| !v.is_empty()).collect();
                coefficients.push(sign * mon[0].parse::<i32>().unwrap_or(1));
                let mut exp: Vec<i32> = vec![0; variables.len()];
                for m in mon.iter() {
                    let tmp: Vec<&str> = m.split('^').collect();
                    let idx = variables.iter().position(|n| n == tmp[0]).unwrap_or(usize::MAX);
                    if idx != usize::MAX {
                        exp[idx] = if tmp.len() == 2 {
                            tmp[1].parse::<i32>().unwrap()
                        } else {
                            1
                        };
                    }
                }
                exponents.append(&mut exp);
                sign = -1;
            }
        }
    }
    return (variables, characteristic, lengths, coefficients, exponents);
}

#[cfg(test)]
mod tests {
    #[test]
    #[should_panic]
    fn test_read_file() {
        let config = super::Config { file_path: "".to_string() };
        super::read_file(config);
    }
}
