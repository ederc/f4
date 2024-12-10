use is_prime::is_prime;

// Calculates (n^x) % p
fn modular_exponent(mut n:i32 ,mut x:i32 , p:i32) -> i32{
    let mut ans = 1;
    if x <= 0 {
        return 1;
    }
    loop {
        if x == 1 {
            return (ans * n) % p;
        }
        if x & 1 == 0 {
            n = (n * n) % p;
            x >>= 1;
            continue;
        } else {
            ans = (ans * n) % p;
            x -= 1;
        }
    }
}

fn modular_inverse(n:i32, p:i32) -> i32{
    debug_assert!(is_prime(&p.to_string()));

    // Return Modular Multiplicative Inverse, that is (n^(p-2)) mod p
    // From Fermat's little theorem
    return modular_exponent(n, p-2, p);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_modular_exponent() {
        assert_eq!(modular_exponent(-2, 5, 7), -4);
        assert_eq!(modular_exponent(-2, 3, 7), -1);
        assert_eq!(modular_exponent(2, 1, 3), modular_inverse(2, 3));
    }

    #[test]
    fn test_modular_inverse() {
        assert_eq!(modular_inverse(2, 3), 2);
        assert_eq!(modular_inverse(1, 2147483647), 1);
        assert_eq!(modular_inverse(-2, 7), -4);
    }
}
