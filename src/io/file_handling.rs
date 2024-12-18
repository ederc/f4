use std::env;
use std::fs;
use is_prime::is_prime;
use crate::types::*;

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
pub fn read_file(config: Config) ->
        (Vec<String>, i32, Vec<usize>, Vec<Vec<i32>>, Vec<Vec<Vec<Exponent>>>) {
    let contents  = fs::read_to_string(config.file_path)
        .expect("Should have been able to read the file")
        .replace(" ", "");
    return read_input_system(contents);
}

fn read_input_system(contents: String) ->
        (Vec<String>, i32, Vec<usize>, Vec<Vec<i32>>, Vec<Vec<Vec<u32>>>) {
    let lines: Vec<&str> = contents.lines().collect();
    let variables: Vec<String> = lines[0].split(",").map(String::from).collect();
    let characteristic: i32 = lines[1].parse::<i32>().unwrap_or(-1);
    if characteristic != 0 && !is_prime(&characteristic.to_string()) {
        panic!("Wrong characteristic!");
    }
    // get correct polynomials, i.e. several polys per line, one poly over several lines, etc.
    let tmp = &lines[2..lines.len()].join("");
    let polynomials: Vec<&str> = tmp.split(",").collect();
    // read in each equation
    let mut exponents: Vec<Vec<Vec<Exponent>>> = Vec::new();
    let mut coefficients: Vec<Vec<i32>>   = Vec::new();
    let mut lengths: Vec<usize>           = Vec::new();
    for p in polynomials {
        let mut cfs: Vec<i32>       = Vec::new();
        let mut exps: Vec<Vec<Exponent>> = Vec::new();
        lengths.push(p.chars().filter(|v| *v == '+' || *v == '-').count() + 1);
        for positive_start in p.split('+') {
            let monomials = positive_start.split('-').filter(|v| !v.is_empty());
            let mut sign: i32 = 1;
            for mons in monomials {
                // first entry is the only positive one
                let mon: Vec<&str> = mons.split('*').filter(|v| !v.is_empty()).collect();
                cfs.push(sign * mon[0].parse::<i32>().unwrap_or(1));
                let mut exp: Vec<Exponent> = vec![0; variables.len()];
                for m in mon.iter() {
                    let tmp: Vec<&str> = m.split('^').collect();
                    let idx = variables.iter().position(|n| n == tmp[0]).unwrap_or(usize::MAX);
                    if idx != usize::MAX {
                        // Applying `+=` and not `=` to also read in monomials
                        // given e.g. as x^2*y*x correctly.
                        exp[idx] += if tmp.len() == 2 {
                            tmp[1].parse::<Exponent>().unwrap()
                        } else {
                            1
                        };
                    }
                }
                exps.push(exp);
                sign = -1;
            }
        }
        exponents.push(exps);
        coefficients.push(cfs);
    }
    return (variables, characteristic, lengths, coefficients, exponents);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn test_read_file_no_input() {
        let config = Config { file_path: "".to_string() };
        read_file(config);
    }

    #[test]
    fn test_read_input_system() {
        let contents = "x,y\n32003\n12*x^2*y-345,100*x-\n3*y^2";
        let (variables, characteristic, lengths, coefficients, exponents)
            = read_input_system(contents.to_string());
        assert_eq!(variables, ["x", "y"]);
        assert_eq!(characteristic, 32003);
        assert_eq!(lengths, [2, 2]);
        assert_eq!(coefficients, [[12, -345], [100, -3]]);
        assert_eq!(exponents, [[[2, 1], [0, 0]], [[1, 0], [0, 2]]]);
    }
    #[test]
    fn test_read_monomials_with_doubled_variables() {
        let contents = "x,y\n32003\n12*x^2*y*y^3";
        let (variables, characteristic, lengths, coefficients, exponents)
            = read_input_system(contents.to_string());
        assert_eq!(variables, ["x", "y"]);
        assert_eq!(characteristic, 32003);
        assert_eq!(lengths, [1]);
        assert_eq!(coefficients, [[12]]);
        assert_eq!(exponents, [[[2, 4]]]);
    }
}
