/// AST for selection expressions.
#[derive(Debug, Clone)]
pub enum Selector {
    /// All atoms.
    All,
    /// No atoms.
    None,
    /// Chain identifier, e.g. `chain A`.
    Chain(char),
    /// Single residue number or range, e.g. `resi 10` or `resi 10-50`.
    Resi(i32, i32),
    /// Atom name, e.g. `name CA`.
    Name(String),
    /// Residue name, e.g. `resn ALA`.
    Resn(String),
    /// Element symbol, e.g. `elem C`.
    Elem(String),
    /// HETATM atoms.
    Hetatm,
    /// Boolean AND.
    And(Box<Selector>, Box<Selector>),
    /// Boolean OR.
    Or(Box<Selector>, Box<Selector>),
    /// Boolean NOT.
    Not(Box<Selector>),
}

/// Parse a selection expression string into a Selector AST.
/// Returns Err with a message if parsing fails.
pub fn parse_selection(input: &str) -> Result<Selector, String> {
    let tokens = tokenize(input)?;
    if tokens.is_empty() {
        return Ok(Selector::All);
    }
    let mut pos = 0;
    let result = parse_or(&tokens, &mut pos)?;
    if pos < tokens.len() {
        return Err(format!("Unexpected token: '{}'", tokens[pos]));
    }
    Ok(result)
}

// ---------- Tokenizer ----------

fn tokenize(input: &str) -> Result<Vec<String>, String> {
    let mut tokens = Vec::new();
    let mut chars = input.chars().peekable();

    while let Some(&c) = chars.peek() {
        if c.is_whitespace() {
            chars.next();
            continue;
        }
        if c == '(' || c == ')' {
            tokens.push(c.to_string());
            chars.next();
            continue;
        }
        // Collect a word (including digits, hyphens for ranges like 10-50)
        let mut word = String::new();
        while let Some(&ch) = chars.peek() {
            if ch.is_whitespace() || ch == '(' || ch == ')' {
                break;
            }
            word.push(ch);
            chars.next();
        }
        if !word.is_empty() {
            tokens.push(word);
        }
    }

    Ok(tokens)
}

// ---------- Recursive descent parser ----------
// Grammar:
//   or_expr   = and_expr ("or" and_expr)*
//   and_expr  = not_expr ("and" not_expr)*
//   not_expr  = "not" not_expr | primary
//   primary   = "(" or_expr ")" | keyword_selector

fn parse_or(tokens: &[String], pos: &mut usize) -> Result<Selector, String> {
    let mut left = parse_and(tokens, pos)?;
    while *pos < tokens.len() && tokens[*pos].eq_ignore_ascii_case("or") {
        *pos += 1;
        let right = parse_and(tokens, pos)?;
        left = Selector::Or(Box::new(left), Box::new(right));
    }
    Ok(left)
}

fn parse_and(tokens: &[String], pos: &mut usize) -> Result<Selector, String> {
    let mut left = parse_not(tokens, pos)?;
    while *pos < tokens.len() && tokens[*pos].eq_ignore_ascii_case("and") {
        *pos += 1;
        let right = parse_not(tokens, pos)?;
        left = Selector::And(Box::new(left), Box::new(right));
    }
    Ok(left)
}

fn parse_not(tokens: &[String], pos: &mut usize) -> Result<Selector, String> {
    if *pos < tokens.len() && tokens[*pos].eq_ignore_ascii_case("not") {
        *pos += 1;
        let inner = parse_not(tokens, pos)?;
        return Ok(Selector::Not(Box::new(inner)));
    }
    parse_primary(tokens, pos)
}

fn parse_primary(tokens: &[String], pos: &mut usize) -> Result<Selector, String> {
    if *pos >= tokens.len() {
        return Err("Unexpected end of selection expression".into());
    }

    let tok = &tokens[*pos];

    // Parenthesized subexpression
    if tok == "(" {
        *pos += 1;
        let inner = parse_or(tokens, pos)?;
        if *pos >= tokens.len() || tokens[*pos] != ")" {
            return Err("Missing closing parenthesis".into());
        }
        *pos += 1;
        return Ok(inner);
    }

    let lower = tok.to_ascii_lowercase();

    match lower.as_str() {
        "all" => {
            *pos += 1;
            Ok(Selector::All)
        }
        "none" => {
            *pos += 1;
            Ok(Selector::None)
        }
        "hetatm" => {
            *pos += 1;
            Ok(Selector::Hetatm)
        }
        "chain" | "c." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "chain")?;
            let ch = arg.chars().next().unwrap_or(' ');
            Ok(Selector::Chain(ch))
        }
        "resi" | "i." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "resi")?;
            parse_resi_range(&arg)
        }
        "name" | "n." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "name")?;
            Ok(Selector::Name(arg.to_uppercase()))
        }
        "resn" | "r." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "resn")?;
            Ok(Selector::Resn(arg.to_uppercase()))
        }
        "elem" | "e." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "elem")?;
            Ok(Selector::Elem(arg))
        }
        _ => {
            Err(format!("Unknown selector keyword: '{}'", tok))
        }
    }
}

fn next_arg(tokens: &[String], pos: &mut usize, keyword: &str) -> Result<String, String> {
    if *pos >= tokens.len() {
        return Err(format!("'{}' requires an argument", keyword));
    }
    let val = tokens[*pos].clone();
    *pos += 1;
    Ok(val)
}

fn parse_resi_range(s: &str) -> Result<Selector, String> {
    if let Some(idx) = s.find('-') {
        // Could be negative number or range
        if idx == 0 {
            // Negative number like -5 or negative range like -5-10
            if let Some(idx2) = s[1..].find('-') {
                let start: i32 = s[..idx2 + 1]
                    .parse()
                    .map_err(|_| format!("Invalid resi range: '{}'", s))?;
                let end: i32 = s[idx2 + 2..]
                    .parse()
                    .map_err(|_| format!("Invalid resi range: '{}'", s))?;
                Ok(Selector::Resi(start, end))
            } else {
                // Just a negative number
                let n: i32 = s.parse().map_err(|_| format!("Invalid resi: '{}'", s))?;
                Ok(Selector::Resi(n, n))
            }
        } else {
            // Normal range like 10-50
            let start: i32 = s[..idx]
                .parse()
                .map_err(|_| format!("Invalid resi range: '{}'", s))?;
            let end: i32 = s[idx + 1..]
                .parse()
                .map_err(|_| format!("Invalid resi range: '{}'", s))?;
            Ok(Selector::Resi(start, end))
        }
    } else {
        let n: i32 = s.parse().map_err(|_| format!("Invalid resi: '{}'", s))?;
        Ok(Selector::Resi(n, n))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all() {
        let sel = parse_selection("all").unwrap();
        assert!(matches!(sel, Selector::All));
    }

    #[test]
    fn test_chain() {
        let sel = parse_selection("chain A").unwrap();
        assert!(matches!(sel, Selector::Chain('A')));
    }

    #[test]
    fn test_resi_range() {
        let sel = parse_selection("resi 10-50").unwrap();
        assert!(matches!(sel, Selector::Resi(10, 50)));
    }

    #[test]
    fn test_boolean() {
        let sel = parse_selection("chain A and resi 1-10").unwrap();
        assert!(matches!(sel, Selector::And(_, _)));
    }

    #[test]
    fn test_not() {
        let sel = parse_selection("not hetatm").unwrap();
        assert!(matches!(sel, Selector::Not(_)));
    }

    #[test]
    fn test_parens() {
        let sel = parse_selection("(chain A or chain B) and name CA").unwrap();
        assert!(matches!(sel, Selector::And(_, _)));
    }

    #[test]
    fn test_empty_is_all() {
        let sel = parse_selection("").unwrap();
        assert!(matches!(sel, Selector::All));
    }
}
