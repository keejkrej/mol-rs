/// AST for selection expressions.
#[derive(Debug, Clone)]
pub enum Selector {
    /// All atoms.
    All,
    /// No atoms.
    None,
    /// Atoms in enabled objects.
    Enabled,
    /// Atoms currently visible through object visibility and representation flags.
    Visible,
    /// Atoms present in the evaluated coordinate state.
    Present,
    /// Atoms with at least one explicit bond.
    Bonded,
    /// Chain identifier, e.g. `chain A`.
    Chain(char),
    /// Chain identifier list or wildcard pattern, e.g. `chain A+B` or `chain A*`.
    ChainPattern(String),
    /// Single residue number or range, e.g. `resi 10` or `resi 10-50`.
    Resi(i32, i32),
    /// Residue numbers and ranges, e.g. `resi 2+12+9+4` or `resi 1+3-5`.
    ResiList(Vec<(i32, i32)>),
    /// Atom name, e.g. `name CA`.
    Name(String),
    /// Residue name, e.g. `resn ALA`.
    Resn(String),
    /// Element symbol, e.g. `elem C`.
    Elem(String),
    /// Alternate location indicator, e.g. `alt A`.
    Alt(char),
    /// Alternate location list or wildcard pattern, e.g. `alt A+B`.
    AltPattern(String),
    /// Secondary structure, e.g. `ss H`, `ss S`, or `ss L`.
    SS(String),
    /// Serial number / atom index, e.g. `serial 10` or `index 10`.
    Serial(i32, i32),
    /// Serial numbers / atom indexes and ranges, e.g. `index 10+20-22`.
    SerialList(Vec<(i32, i32)>),
    /// Object/model name, e.g. `model obj01` or `object ligand*`.
    Model(String),
    /// Visible representation name, e.g. `rep lines` or `rep sticks+spheres`.
    Rep(String),
    /// Atom display color, e.g. `color red`.
    Color(String),
    /// Atom or coordinate property comparison, e.g. `b < 20`, `q >= 0.5`, or `x > 10`.
    Property(AtomProperty, CompareOp, f32),
    /// HETATM atoms.
    Hetatm,
    /// Hydrogen atoms.
    Hydrogen,
    /// Solvent atoms.
    Solvent,
    /// Protein and nucleic polymer atoms.
    Polymer,
    /// Protein polymer atoms.
    Protein,
    /// Nucleic polymer atoms.
    Nucleic,
    /// Metal atoms.
    Metals,
    /// Polymer guide atoms, e.g. CA for protein and C4'/C4* for nucleic residues.
    Guide,
    /// Organic molecules (heuristic: non-HETATM).
    Organic,
    /// Inorganic molecules (heuristic: HETATM).
    Inorganic,
    /// Boolean AND.
    And(Box<Selector>, Box<Selector>),
    /// Boolean OR.
    Or(Box<Selector>, Box<Selector>),
    /// Boolean NOT.
    Not(Box<Selector>),
    /// Select all atoms in residues containing atoms matching inner selection.
    Byres(Box<Selector>),
    /// Select all atoms in chains containing atoms matching inner selection.
    Bychain(Box<Selector>),
    /// Select all atoms in objects containing atoms matching inner selection.
    Byobject(Box<Selector>),
    /// Select all bonded-connected atoms containing atoms matching inner selection.
    Bymolecule(Box<Selector>),
    /// First atom from another selection.
    First(Box<Selector>),
    /// Last atom from another selection.
    Last(Box<Selector>),
    /// Atoms directly bonded to atoms in another selection.
    Neighbor(Box<Selector>),
    /// Atoms around another selection, excluding the inner selection.
    Around(f32, Box<Selector>),
    /// Atoms within distance of another selection.
    Within(f32, Box<Selector>),
    /// Atoms within distance of another selection, including the inner selection.
    Expand(f32, Box<Selector>),
    /// Atoms within a bond-count extension of another selection.
    Extend(usize, Box<Selector>),
    /// Atoms outside distance of another selection.
    Beyond(f32, Box<Selector>),
    /// Atoms within distance of another selection, excluding the inner selection.
    NearTo(f32, Box<Selector>),
}

#[derive(Debug, Clone, Copy)]
pub enum AtomProperty {
    BFactor,
    Occupancy,
    X,
    Y,
    Z,
}

#[derive(Debug, Clone, Copy)]
pub enum CompareOp {
    Greater,
    Less,
    Equal,
    GreaterEqual,
    LessEqual,
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
        if matches!(c, '(' | ')' | '&' | '|' | '!') {
            tokens.push(c.to_string());
            chars.next();
            continue;
        }
        if matches!(c, '>' | '<' | '=') {
            let mut op = String::new();
            op.push(c);
            chars.next();
            if chars.peek() == Some(&'=') {
                op.push('=');
                chars.next();
            }
            tokens.push(op);
            continue;
        }
        // Collect a word (including digits, hyphens for ranges like 10-50).
        let mut word = String::new();
        while let Some(&ch) = chars.peek() {
            if ch.is_whitespace() || matches!(ch, '(' | ')' | '&' | '|' | '!' | '>' | '<' | '=') {
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
    while *pos < tokens.len() && is_or_operator(&tokens[*pos]) {
        *pos += 1;
        let right = parse_and(tokens, pos)?;
        left = Selector::Or(Box::new(left), Box::new(right));
    }
    Ok(left)
}

fn parse_and(tokens: &[String], pos: &mut usize) -> Result<Selector, String> {
    let mut left = parse_not(tokens, pos)?;
    while *pos < tokens.len() && is_and_operator(&tokens[*pos]) {
        *pos += 1;
        let right = parse_not(tokens, pos)?;
        left = Selector::And(Box::new(left), Box::new(right));
    }
    Ok(left)
}

fn parse_not(tokens: &[String], pos: &mut usize) -> Result<Selector, String> {
    if *pos < tokens.len() && is_not_operator(&tokens[*pos]) {
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
        "all" | "*" => {
            *pos += 1;
            Ok(Selector::All)
        }
        "none" => {
            *pos += 1;
            Ok(Selector::None)
        }
        "enabled" => {
            *pos += 1;
            Ok(Selector::Enabled)
        }
        "visible" | "vis" | "v;" | "v." => {
            *pos += 1;
            Ok(Selector::Visible)
        }
        "present" | "pr." => {
            *pos += 1;
            Ok(Selector::Present)
        }
        "bonded" => {
            *pos += 1;
            Ok(Selector::Bonded)
        }
        "hetatm" | "het" => {
            *pos += 1;
            Ok(Selector::Hetatm)
        }
        "hydrogens" | "hydro" | "h;" | "h." => {
            *pos += 1;
            Ok(Selector::Hydrogen)
        }
        "solvent" | "sol." => {
            *pos += 1;
            Ok(Selector::Solvent)
        }
        "polymer" | "pol." => {
            *pos += 1;
            Ok(Selector::Polymer)
        }
        "polymer.protein" | "protein" | "pro." => {
            *pos += 1;
            Ok(Selector::Protein)
        }
        "polymer.nucleic" | "nucleic" | "nuc." => {
            *pos += 1;
            Ok(Selector::Nucleic)
        }
        "metals" => {
            *pos += 1;
            Ok(Selector::Metals)
        }
        "guide" => {
            *pos += 1;
            Ok(Selector::Guide)
        }
        "organic" | "org." => {
            *pos += 1;
            Ok(Selector::Organic)
        }
        "inorganic" | "ino." => {
            *pos += 1;
            Ok(Selector::Inorganic)
        }
        "chain" | "c;" | "c." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "chain")?;
            let mut chars = arg.chars();
            if let (Some(ch), None) = (chars.next(), chars.next()) {
                Ok(Selector::Chain(ch))
            } else {
                Ok(Selector::ChainPattern(arg))
            }
        }
        "resi" | "i;" | "i." | "residue" | "resident" | "resid" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "resi")?;
            parse_resi_range(&arg)
        }
        "name" | "n;" | "n." | "ca" => {
            *pos += 1;
            if lower.as_str() == "ca" {
                Ok(Selector::Name("CA".to_string()))
            } else {
                let arg = next_arg(tokens, pos, "name")?;
                Ok(Selector::Name(arg.to_uppercase()))
            }
        }
        "resn" | "r;" | "r." | "resname" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "resn")?;
            Ok(Selector::Resn(arg.to_uppercase()))
        }
        "serial" | "index" | "idx." | "id" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "serial")?;
            parse_serial_range(&arg)
        }
        "object" | "model" | "o." | "m;" | "m." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "model")?;
            Ok(Selector::Model(arg))
        }
        "rep" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "rep")?;
            Ok(Selector::Rep(arg))
        }
        "color" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "color")?;
            Ok(Selector::Color(arg))
        }
        "b" => {
            *pos += 1;
            parse_property_selector(AtomProperty::BFactor, tokens, pos, "b")
        }
        "q" => {
            *pos += 1;
            parse_property_selector(AtomProperty::Occupancy, tokens, pos, "q")
        }
        "x" => {
            *pos += 1;
            parse_property_selector(AtomProperty::X, tokens, pos, "x")
        }
        "y" => {
            *pos += 1;
            parse_property_selector(AtomProperty::Y, tokens, pos, "y")
        }
        "z" => {
            *pos += 1;
            parse_property_selector(AtomProperty::Z, tokens, pos, "z")
        }
        "elem" | "e;" | "e." | "element" | "symbol" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "elem")?;
            Ok(Selector::Elem(arg))
        }
        "alt" | "altloc" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "alt")?;
            if let Some(ch) = parse_special_char_arg(&arg) {
                Ok(Selector::Alt(ch))
            } else {
                let mut chars = arg.chars();
                if let (Some(ch), None) = (chars.next(), chars.next()) {
                    Ok(Selector::Alt(ch))
                } else {
                    Ok(Selector::AltPattern(arg))
                }
            }
        }
        "ss" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "ss")?;
            Ok(Selector::SS(arg.to_uppercase()))
        }
        "within" | "w." => {
            *pos += 1;
            let distance = parse_distance_arg(tokens, pos, "within")?;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Within(distance, Box::new(inner)))
        }
        "around" | "a;" | "a." => {
            *pos += 1;
            let distance = parse_distance_arg(tokens, pos, "around")?;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Around(distance, Box::new(inner)))
        }
        "expand" | "x;" | "x." => {
            *pos += 1;
            let distance = parse_distance_arg(tokens, pos, "expand")?;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Expand(distance, Box::new(inner)))
        }
        "extend" | "xt." => {
            *pos += 1;
            let count_token = next_arg(tokens, pos, "extend")?;
            let count: usize = count_token
                .parse()
                .map_err(|_| format!("Invalid bond count '{}'", count_token))?;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Extend(count, Box::new(inner)))
        }
        "beyond" | "be." => {
            *pos += 1;
            let distance = parse_distance_arg(tokens, pos, "beyond")?;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Beyond(distance, Box::new(inner)))
        }
        "near_to" | "nto." => {
            *pos += 1;
            let distance = parse_distance_arg(tokens, pos, "near_to")?;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::NearTo(distance, Box::new(inner)))
        }
        "byresidue" | "byresi" | "byres" | "br;" | "br." => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Byres(Box::new(inner)))
        }
        "byca" | "bycalpha" | "bca." => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::And(
                Box::new(Selector::Byres(Box::new(inner))),
                Box::new(Selector::Name("CA".to_string())),
            ))
        }
        "bychain" | "bc." => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Bychain(Box::new(inner)))
        }
        "byobject" | "byobj" | "bo;" | "bo." => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Byobject(Box::new(inner)))
        }
        "bymolecule" | "bymol" | "bm." => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Bymolecule(Box::new(inner)))
        }
        "first" => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::First(Box::new(inner)))
        }
        "last" => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Last(Box::new(inner)))
        }
        "neighbor" | "nbr" => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Neighbor(Box::new(inner)))
        }
        "bb" | "bb." | "backbone" => {
            *pos += 1;
            Ok(Selector::Or(
                Box::new(Selector::Or(
                    Box::new(Selector::Name("N".to_string())),
                    Box::new(Selector::Name("CA".to_string())),
                )),
                Box::new(Selector::Or(
                    Box::new(Selector::Name("C".to_string())),
                    Box::new(Selector::Name("O".to_string())),
                )),
            ))
        }
        "sc" | "sc." | "sidechain" => {
            *pos += 1;
            Ok(Selector::Not(Box::new(backbone_selector())))
        }
        _ => Err(format!("Unknown selector keyword: '{}'", tok)),
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

fn is_and_operator(token: &str) -> bool {
    token == "&" || token.eq_ignore_ascii_case("and")
}

fn is_or_operator(token: &str) -> bool {
    token == "|" || token.eq_ignore_ascii_case("or")
}

fn is_not_operator(token: &str) -> bool {
    token == "!" || token.eq_ignore_ascii_case("not")
}

fn parse_property_selector(
    property: AtomProperty,
    tokens: &[String],
    pos: &mut usize,
    keyword: &str,
) -> Result<Selector, String> {
    let op_token = next_arg(tokens, pos, keyword)?;
    let value_token = next_arg(tokens, pos, keyword)?;
    let op = parse_compare_op(&op_token)
        .ok_or_else(|| format!("Invalid comparison operator for '{keyword}': '{op_token}'"))?;
    let value: f32 = value_token
        .parse()
        .map_err(|_| format!("Invalid comparison value for '{keyword}': '{value_token}'"))?;

    Ok(Selector::Property(property, op, value))
}

fn parse_distance_arg(tokens: &[String], pos: &mut usize, keyword: &str) -> Result<f32, String> {
    let distance_token = next_arg(tokens, pos, keyword)?;
    distance_token
        .parse()
        .map_err(|_| format!("Invalid distance '{}'", distance_token))
}

fn parse_unary_inner(tokens: &[String], pos: &mut usize) -> Result<Selector, String> {
    if *pos < tokens.len() && tokens[*pos] == "(" {
        parse_primary(tokens, pos)
    } else {
        parse_not(tokens, pos)
    }
}

fn parse_compare_op(s: &str) -> Option<CompareOp> {
    match s {
        ">" => Some(CompareOp::Greater),
        "<" => Some(CompareOp::Less),
        "=" | "==" => Some(CompareOp::Equal),
        ">=" => Some(CompareOp::GreaterEqual),
        "<=" => Some(CompareOp::LessEqual),
        _ => None,
    }
}

fn parse_special_char_arg(s: &str) -> Option<char> {
    match s {
        "\"\"" | "''" | "blank" | "none" => Some(' '),
        _ => None,
    }
}

fn backbone_selector() -> Selector {
    Selector::Or(
        Box::new(Selector::Or(
            Box::new(Selector::Name("N".to_string())),
            Box::new(Selector::Name("CA".to_string())),
        )),
        Box::new(Selector::Or(
            Box::new(Selector::Name("C".to_string())),
            Box::new(Selector::Name("O".to_string())),
        )),
    )
}

fn parse_resi_range(s: &str) -> Result<Selector, String> {
    if s.contains('+') {
        Ok(Selector::ResiList(parse_numeric_range_list(s, "resi")?))
    } else {
        let (start, end) = parse_numeric_range(s, "resi")?;
        Ok(Selector::Resi(start, end))
    }
}

fn parse_serial_range(s: &str) -> Result<Selector, String> {
    if s.contains('+') {
        Ok(Selector::SerialList(parse_numeric_range_list(s, "serial")?))
    } else {
        let (start, end) = parse_numeric_range(s, "serial")?;
        Ok(Selector::Serial(start, end))
    }
}

fn parse_numeric_range_list(s: &str, label: &str) -> Result<Vec<(i32, i32)>, String> {
    s.split('+')
        .map(|part| parse_numeric_range(part.trim(), label))
        .collect()
}

fn parse_numeric_range(s: &str, label: &str) -> Result<(i32, i32), String> {
    if let Some(idx) = s.find('-') {
        if idx == 0 {
            if let Some(idx2) = s[1..].find('-') {
                let start: i32 = s[..idx2 + 1]
                    .parse()
                    .map_err(|_| format!("Invalid {label} range: '{s}'"))?;
                let end: i32 = s[idx2 + 2..]
                    .parse()
                    .map_err(|_| format!("Invalid {label} range: '{s}'"))?;
                Ok((start, end))
            } else {
                let n: i32 = s.parse().map_err(|_| format!("Invalid {label}: '{s}'"))?;
                Ok((n, n))
            }
        } else {
            let start: i32 = s[..idx]
                .parse()
                .map_err(|_| format!("Invalid {label} range: '{s}'"))?;
            let end: i32 = s[idx + 1..]
                .parse()
                .map_err(|_| format!("Invalid {label} range: '{s}'"))?;
            Ok((start, end))
        }
    } else {
        let n: i32 = s.parse().map_err(|_| format!("Invalid {label}: '{s}'"))?;
        Ok((n, n))
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
    fn test_chain_pattern_selection() {
        assert!(matches!(
            parse_selection("chain A+B").unwrap(),
            Selector::ChainPattern(ref s) if s == "A+B"
        ));
        assert!(matches!(
            parse_selection("chain A*").unwrap(),
            Selector::ChainPattern(ref s) if s == "A*"
        ));
    }

    #[test]
    fn test_resi_range() {
        let sel = parse_selection("resi 10-50").unwrap();
        assert!(matches!(sel, Selector::Resi(10, 50)));
    }

    #[test]
    fn test_plus_separated_alpha_list_args() {
        assert!(matches!(
            parse_selection("name C+N+O+CA").unwrap(),
            Selector::Name(ref s) if s == "C+N+O+CA"
        ));
        assert!(matches!(
            parse_selection("resn ALA+GLY").unwrap(),
            Selector::Resn(ref s) if s == "ALA+GLY"
        ));
        assert!(matches!(
            parse_selection("elem C+N").unwrap(),
            Selector::Elem(ref s) if s == "C+N"
        ));
    }

    #[test]
    fn test_plus_separated_numeric_list_args() {
        assert!(matches!(
            parse_selection("resi 2+12+9+4").unwrap(),
            Selector::ResiList(ref ranges)
                if ranges == &vec![(2, 2), (12, 12), (9, 9), (4, 4)]
        ));
        assert!(matches!(
            parse_selection("resi 1+3-5").unwrap(),
            Selector::ResiList(ref ranges) if ranges == &vec![(1, 1), (3, 5)]
        ));
        assert!(matches!(
            parse_selection("index 10+20-22").unwrap(),
            Selector::SerialList(ref ranges) if ranges == &vec![(10, 10), (20, 22)]
        ));
    }

    #[test]
    fn test_selector_aliases() {
        assert!(matches!(parse_selection("*").unwrap(), Selector::All));
        assert!(matches!(
            parse_selection("residue 10").unwrap(),
            Selector::Resi(10, 10)
        ));
        assert!(matches!(
            parse_selection("i; 10").unwrap(),
            Selector::Resi(10, 10)
        ));
        assert!(matches!(
            parse_selection("resname ALA").unwrap(),
            Selector::Resn(ref s) if s == "ALA"
        ));
        assert!(matches!(
            parse_selection("r; ALA").unwrap(),
            Selector::Resn(ref s) if s == "ALA"
        ));
        assert!(matches!(
            parse_selection("element C").unwrap(),
            Selector::Elem(_)
        ));
        assert!(matches!(
            parse_selection("e; C").unwrap(),
            Selector::Elem(_)
        ));
        assert!(matches!(
            parse_selection("symbol C").unwrap(),
            Selector::Elem(_)
        ));
        assert!(matches!(
            parse_selection("alt A").unwrap(),
            Selector::Alt('A')
        ));
        assert!(matches!(
            parse_selection("alt A+B").unwrap(),
            Selector::AltPattern(ref s) if s == "A+B"
        ));
        assert!(matches!(
            parse_selection("alt blank").unwrap(),
            Selector::Alt(' ')
        ));
        assert!(matches!(
            parse_selection("ss H").unwrap(),
            Selector::SS(ref s) if s == "H"
        ));
        assert!(matches!(
            parse_selection("ss H+S").unwrap(),
            Selector::SS(ref s) if s == "H+S"
        ));
        assert!(matches!(
            parse_selection("ca").unwrap(),
            Selector::Name(ref s) if s == "CA"
        ));
        assert!(matches!(
            parse_selection("n; CA").unwrap(),
            Selector::Name(ref s) if s == "CA"
        ));
        assert!(matches!(parse_selection("het").unwrap(), Selector::Hetatm));
        assert!(matches!(
            parse_selection("hydrogens").unwrap(),
            Selector::Hydrogen
        ));
        assert!(matches!(
            parse_selection("hydro").unwrap(),
            Selector::Hydrogen
        ));
        assert!(matches!(parse_selection("h.").unwrap(), Selector::Hydrogen));
        assert!(matches!(parse_selection("h;").unwrap(), Selector::Hydrogen));
        assert!(matches!(
            parse_selection("solvent").unwrap(),
            Selector::Solvent
        ));
        assert!(matches!(
            parse_selection("sol.").unwrap(),
            Selector::Solvent
        ));
        assert!(matches!(
            parse_selection("polymer").unwrap(),
            Selector::Polymer
        ));
        assert!(matches!(
            parse_selection("pol.").unwrap(),
            Selector::Polymer
        ));
        assert!(matches!(
            parse_selection("polymer.protein").unwrap(),
            Selector::Protein
        ));
        assert!(matches!(
            parse_selection("protein").unwrap(),
            Selector::Protein
        ));
        assert!(matches!(
            parse_selection("pro.").unwrap(),
            Selector::Protein
        ));
        assert!(matches!(
            parse_selection("polymer.nucleic").unwrap(),
            Selector::Nucleic
        ));
        assert!(matches!(
            parse_selection("nucleic").unwrap(),
            Selector::Nucleic
        ));
        assert!(matches!(
            parse_selection("nuc.").unwrap(),
            Selector::Nucleic
        ));
        assert!(matches!(
            parse_selection("metals").unwrap(),
            Selector::Metals
        ));
        assert!(matches!(
            parse_selection("organic").unwrap(),
            Selector::Organic
        ));
        assert!(matches!(
            parse_selection("org.").unwrap(),
            Selector::Organic
        ));
        assert!(matches!(
            parse_selection("inorganic").unwrap(),
            Selector::Inorganic
        ));
        assert!(matches!(
            parse_selection("ino.").unwrap(),
            Selector::Inorganic
        ));
        assert!(matches!(
            parse_selection("enabled").unwrap(),
            Selector::Enabled
        ));
        assert!(matches!(
            parse_selection("visible").unwrap(),
            Selector::Visible
        ));
        assert!(matches!(parse_selection("vis").unwrap(), Selector::Visible));
        assert!(matches!(parse_selection("v;").unwrap(), Selector::Visible));
        assert!(matches!(parse_selection("v.").unwrap(), Selector::Visible));
        assert!(matches!(
            parse_selection("present").unwrap(),
            Selector::Present
        ));
        assert!(matches!(parse_selection("pr.").unwrap(), Selector::Present));
        assert!(matches!(
            parse_selection("bonded").unwrap(),
            Selector::Bonded
        ));
        assert!(matches!(parse_selection("guide").unwrap(), Selector::Guide));
    }

    #[test]
    fn test_serial_selection() {
        assert!(matches!(
            parse_selection("serial 5").unwrap(),
            Selector::Serial(5, 5)
        ));
        assert!(matches!(
            parse_selection("index 10-20").unwrap(),
            Selector::Serial(10, 20)
        ));
        assert!(matches!(
            parse_selection("idx. 10-20").unwrap(),
            Selector::Serial(10, 20)
        ));
    }

    #[test]
    fn test_model_selection() {
        assert!(matches!(
            parse_selection("model prot*").unwrap(),
            Selector::Model(ref s) if s == "prot*"
        ));
        assert!(matches!(
            parse_selection("object ligand").unwrap(),
            Selector::Model(ref s) if s == "ligand"
        ));
        assert!(matches!(
            parse_selection("m. ligand").unwrap(),
            Selector::Model(ref s) if s == "ligand"
        ));
        assert!(matches!(
            parse_selection("m; ligand").unwrap(),
            Selector::Model(ref s) if s == "ligand"
        ));
        assert!(matches!(
            parse_selection("o. ligand").unwrap(),
            Selector::Model(ref s) if s == "ligand"
        ));
    }

    #[test]
    fn test_rep_selection() {
        assert!(matches!(
            parse_selection("rep lines").unwrap(),
            Selector::Rep(ref s) if s == "lines"
        ));
        assert!(matches!(
            parse_selection("rep sticks+spheres").unwrap(),
            Selector::Rep(ref s) if s == "sticks+spheres"
        ));
    }

    #[test]
    fn test_color_selection() {
        assert!(matches!(
            parse_selection("color red").unwrap(),
            Selector::Color(ref s) if s == "red"
        ));
        assert!(matches!(
            parse_selection("color grey").unwrap(),
            Selector::Color(ref s) if s == "grey"
        ));
    }

    #[test]
    fn test_property_selection() {
        assert!(matches!(
            parse_selection("b < 20").unwrap(),
            Selector::Property(AtomProperty::BFactor, CompareOp::Less, v)
                if (v - 20.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("q>=0.5").unwrap(),
            Selector::Property(AtomProperty::Occupancy, CompareOp::GreaterEqual, v)
                if (v - 0.5).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("x == -1.25").unwrap(),
            Selector::Property(AtomProperty::X, CompareOp::Equal, v)
                if (v + 1.25).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("z<=10").unwrap(),
            Selector::Property(AtomProperty::Z, CompareOp::LessEqual, v)
                if (v - 10.0).abs() < f32::EPSILON
        ));
    }

    #[test]
    fn test_within_selector() {
        assert!(matches!(
            parse_selection("within 4 chain A").unwrap(),
            Selector::Within(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("w. 4 chain A").unwrap(),
            Selector::Within(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
    }

    #[test]
    fn test_around_expand_extend_selectors() {
        assert!(matches!(
            parse_selection("a. 4 chain A").unwrap(),
            Selector::Around(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("a; 4 chain A").unwrap(),
            Selector::Around(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("expand 4 chain A").unwrap(),
            Selector::Expand(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("x. 4 chain A").unwrap(),
            Selector::Expand(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("x; 4 chain A").unwrap(),
            Selector::Expand(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("extend 2 chain A").unwrap(),
            Selector::Extend(2, _)
        ));
        assert!(matches!(
            parse_selection("xt. 2 chain A").unwrap(),
            Selector::Extend(2, _)
        ));
    }

    #[test]
    fn test_distance_selector_aliases() {
        assert!(matches!(
            parse_selection("beyond 4 chain A").unwrap(),
            Selector::Beyond(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("be. 4 chain A").unwrap(),
            Selector::Beyond(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("near_to 4 chain A").unwrap(),
            Selector::NearTo(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("nto. 4 chain A").unwrap(),
            Selector::NearTo(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
    }

    #[test]
    fn test_within_precedence() {
        let sel = parse_selection("within 4 chain A and name CA").unwrap();
        if let Selector::And(_, right) = sel {
            assert!(matches!(*right, Selector::Name(ref s) if s == "CA"));
        } else {
            panic!("expected 'within' to be grouped before outer 'and'");
        }
    }

    #[test]
    fn test_neighbor_selector() {
        assert!(matches!(
            parse_selection("neighbor chain A").unwrap(),
            Selector::Neighbor(_)
        ));
    }

    #[test]
    fn test_byres_selector() {
        assert!(matches!(
            parse_selection("byres chain A").unwrap(),
            Selector::Byres(_)
        ));
        assert!(matches!(
            parse_selection("byresidue chain A").unwrap(),
            Selector::Byres(_)
        ));
        assert!(matches!(
            parse_selection("byresi chain A").unwrap(),
            Selector::Byres(_)
        ));
        assert!(matches!(
            parse_selection("br. chain A").unwrap(),
            Selector::Byres(_)
        ));
        assert!(matches!(
            parse_selection("br; chain A").unwrap(),
            Selector::Byres(_)
        ));
        assert!(matches!(
            parse_selection("byca chain A").unwrap(),
            Selector::And(_, _)
        ));
        assert!(matches!(
            parse_selection("bycalpha chain A").unwrap(),
            Selector::And(_, _)
        ));
        assert!(matches!(
            parse_selection("bca. chain A").unwrap(),
            Selector::And(_, _)
        ));
    }

    #[test]
    fn test_by_expansion_selectors() {
        assert!(matches!(
            parse_selection("bychain name CA").unwrap(),
            Selector::Bychain(_)
        ));
        assert!(matches!(
            parse_selection("bc. name CA").unwrap(),
            Selector::Bychain(_)
        ));
        assert!(matches!(
            parse_selection("byobject chain A").unwrap(),
            Selector::Byobject(_)
        ));
        assert!(matches!(
            parse_selection("byobj chain A").unwrap(),
            Selector::Byobject(_)
        ));
        assert!(matches!(
            parse_selection("bo. chain A").unwrap(),
            Selector::Byobject(_)
        ));
        assert!(matches!(
            parse_selection("bo; chain A").unwrap(),
            Selector::Byobject(_)
        ));
        assert!(matches!(
            parse_selection("bymolecule serial 10").unwrap(),
            Selector::Bymolecule(_)
        ));
        assert!(matches!(
            parse_selection("bymol serial 10").unwrap(),
            Selector::Bymolecule(_)
        ));
        assert!(matches!(
            parse_selection("bm. serial 10").unwrap(),
            Selector::Bymolecule(_)
        ));
    }

    #[test]
    fn test_first_last_selectors() {
        assert!(matches!(
            parse_selection("first chain A").unwrap(),
            Selector::First(_)
        ));
        assert!(matches!(
            parse_selection("last chain A").unwrap(),
            Selector::Last(_)
        ));
    }

    #[test]
    fn test_bb_selector() {
        assert!(matches!(parse_selection("bb").unwrap(), Selector::Or(_, _)));
        assert!(matches!(
            parse_selection("backbone").unwrap(),
            Selector::Or(_, _)
        ));
        assert!(matches!(
            parse_selection("sidechain").unwrap(),
            Selector::Not(_)
        ));
    }

    #[test]
    fn test_boolean() {
        let sel = parse_selection("chain A and resi 1-10").unwrap();
        assert!(matches!(sel, Selector::And(_, _)));

        let sel = parse_selection("chain A&resi 1-10").unwrap();
        assert!(matches!(sel, Selector::And(_, _)));

        let sel = parse_selection("chain A | chain B").unwrap();
        assert!(matches!(sel, Selector::Or(_, _)));
    }

    #[test]
    fn test_not() {
        let sel = parse_selection("not hetatm").unwrap();
        assert!(matches!(sel, Selector::Not(_)));

        let sel = parse_selection("!hetatm").unwrap();
        assert!(matches!(sel, Selector::Not(_)));

        let sel = parse_selection("chain A &! name C").unwrap();
        if let Selector::And(_, right) = sel {
            assert!(matches!(*right, Selector::Not(_)));
        } else {
            panic!("expected '&!' to parse as AND followed by NOT");
        }
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
