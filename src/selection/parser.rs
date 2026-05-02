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
    /// Atoms in a stored named selection, e.g. `%sele`.
    Named(String),
    /// Atoms in a stored named selection or object, e.g. `sele` or `object1`.
    Identifier(String),
    /// Atoms present in the evaluated coordinate state.
    Present,
    /// Atoms present in a coordinate state; -1 means current evaluated state.
    State(isize),
    /// Atoms with at least one explicit bond.
    Bonded,
    /// Hydrogen bond donor atoms (derived from current atom/bond data).
    Donors,
    /// Hydrogen bond acceptor atoms (derived from current atom/bond data).
    Acceptors,
    /// Atoms participating in delocalized/aromatic bonds.
    Delocalized,
    /// Atoms with a PyMOL atom flag bit set, e.g. `flag 25`.
    Flag(u8),
    /// Atoms masked from picking/selection.
    Masked,
    /// Atoms protected from movement.
    Protected,
    /// Chain identifier, e.g. `chain A`.
    Chain(char),
    /// Chain identifier list or wildcard pattern, e.g. `chain A+B` or `chain A*`.
    ChainPattern(String),
    /// Segment identifier list or wildcard pattern, e.g. `segi A+B` or `segment PRO*`.
    Segi(String),
    /// Single residue number or range, with optional insertion code boundaries, e.g. `resi 10`, `resi 10-50`, or `resi 9A-10A`.
    Resi(i32, i32, Option<char>, Option<char>),
    /// Residue numbers and ranges, e.g. `resi 2+12+9+4`, `resi 1+3-5`, or `resi 9A+10`.
    ResiList(Vec<(i32, i32, Option<char>, Option<char>)>),
    /// Atom name, e.g. `name CA`.
    Name(String),
    /// Residue name, e.g. `resn ALA`.
    Resn(String),
    /// Peptide sequence pattern, e.g. `pepseq AG`, `ps. A+G`, or `pepseq A-G`.
    Pepseq(String),
    /// Force-field text type, e.g. `text_type C3`.
    TextType(String),
    /// Force-field numeric type, e.g. `numeric_type 42` or `nt. 10+20-30`.
    NumericType(Vec<(i32, i32)>),
    /// Custom atom property string, e.g. `custom ligand`.
    Custom(String),
    /// Atom label text, e.g. `label site*`.
    Label(String),
    /// Atom stereochemistry, e.g. `stereo R` or `stereo odd`.
    Stereo(String),
    /// Element symbol, e.g. `elem C`.
    Elem(String),
    /// Alternate location indicator, e.g. `alt A`.
    Alt(char),
    /// Alternate location list or wildcard pattern, e.g. `alt A+B`.
    AltPattern(String),
    /// Secondary structure, e.g. `ss H`, `ss S`, or `ss L`.
    SS(String),
    /// PDB serial / atom ID, e.g. `serial 10` or `id 10`.
    Serial(i32, i32),
    /// PDB serial / atom ID numbers and ranges, e.g. `id 10+20-22`.
    SerialList(Vec<(i32, i32)>),
    /// 1-based object atom index, e.g. `index 10`.
    Index(i32, i32),
    /// 1-based object atom indexes and ranges, e.g. `idx. 10+20-22`.
    IndexList(Vec<(i32, i32)>),
    /// Original atom load rank, e.g. `rank 10`.
    Rank(i32, i32),
    /// Original atom load ranks and ranges, e.g. `rank 10+20-22`.
    RankList(Vec<(i32, i32)>),
    /// Object/model name, e.g. `model obj01` or `object ligand*`.
    Model(String),
    /// Visible representation name, e.g. `rep lines` or `rep sticks+spheres`.
    Rep(String),
    /// Atom display color, e.g. `color red`.
    Color(String),
    /// Explicit per-atom cartoon color override, e.g. `cartoon_color red`.
    CartoonColor(String),
    /// Explicit per-atom ribbon color override, e.g. `ribbon_color blue`.
    RibbonColor(String),
    /// Atom or coordinate property comparison, e.g. `b < 20`, `q >= 0.5`, or `x > 10`.
    Property(AtomProperty, CompareOp, f32),
    /// Named custom atom property comparison, e.g. `p.score > 0.5` or `p.kind in ligand*`.
    CustomProperty(String, CustomPropertyOp, String),
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
    /// PyMOL `left in right`: left atoms whose full atom identity appears in right.
    In(Box<Selector>, Box<Selector>),
    /// PyMOL `left like right`: left atoms whose residue number/name identity appears in right.
    Like(Box<Selector>, Box<Selector>),
    /// Select all atoms in residues containing atoms matching inner selection.
    Byres(Box<Selector>),
    /// Select all atoms in chains containing atoms matching inner selection.
    Bychain(Box<Selector>),
    /// Select all atoms in segments containing atoms matching inner selection.
    Bysegment(Box<Selector>),
    /// Select all atoms in objects containing atoms matching inner selection.
    Byobject(Box<Selector>),
    /// Select all bonded-connected atoms containing atoms matching inner selection.
    Bymolecule(Box<Selector>),
    /// Select all atoms in rings containing atoms matching inner selection.
    Byring(Box<Selector>),
    /// First atom from another selection.
    First(Box<Selector>),
    /// Last atom from another selection.
    Last(Box<Selector>),
    /// Atoms directly bonded to atoms in another selection.
    Neighbor(Box<Selector>),
    /// Atoms directly bonded to atoms in another selection, including selected atoms bonded to selected atoms.
    BoundTo(Box<Selector>),
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
    /// Atoms separated from another selection by at least distance plus VDW radii.
    Gap(f32, Box<Selector>),
}

#[derive(Debug, Clone, Copy)]
pub enum AtomProperty {
    BFactor,
    Occupancy,
    FormalCharge,
    PartialCharge,
    Vdw,
    ElecRadius,
    Cartoon,
    Geom,
    Valence,
    Reps,
    Protons,
    Flags,
    ExplicitDegree,
    ExplicitValence,
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

#[derive(Debug, Clone, Copy)]
pub enum CustomPropertyOp {
    Greater,
    Less,
    Equal,
    GreaterEqual,
    LessEqual,
    In,
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
//   or_expr   = and_expr (("or" | "+" | "in" | "like") and_expr)*
//   and_expr  = not_expr (("and" | "-") not_expr)*
//   not_expr  = "not" not_expr | primary
//   primary   = "(" or_expr ")" | keyword_selector

fn parse_or(tokens: &[String], pos: &mut usize) -> Result<Selector, String> {
    let mut left = parse_and(tokens, pos)?;
    while *pos < tokens.len()
        && (is_or_operator(&tokens[*pos])
            || is_in_operator(&tokens[*pos])
            || is_like_operator(&tokens[*pos]))
    {
        let op = tokens[*pos].clone();
        *pos += 1;
        let right = parse_and(tokens, pos)?;
        left = if is_in_operator(&op) {
            Selector::In(Box::new(left), Box::new(right))
        } else if is_like_operator(&op) {
            Selector::Like(Box::new(left), Box::new(right))
        } else {
            Selector::Or(Box::new(left), Box::new(right))
        };
    }
    Ok(left)
}

fn parse_and(tokens: &[String], pos: &mut usize) -> Result<Selector, String> {
    let mut left = parse_not(tokens, pos)?;
    while *pos < tokens.len()
        && (is_and_operator(&tokens[*pos]) || is_subtract_operator(&tokens[*pos]))
    {
        let subtract = is_subtract_operator(&tokens[*pos]);
        *pos += 1;
        let right = parse_not(tokens, pos)?;
        left = if subtract {
            Selector::And(Box::new(left), Box::new(Selector::Not(Box::new(right))))
        } else {
            Selector::And(Box::new(left), Box::new(right))
        };
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
    if tok == "%" {
        *pos += 1;
        let name = next_arg(tokens, pos, "%")?;
        return Ok(Selector::Named(name));
    }
    if let Some(name) = tok.strip_prefix('%') {
        if !name.is_empty() {
            *pos += 1;
            return Ok(Selector::Named(name.to_string()));
        }
    }
    if lower == "p." {
        *pos += 1;
        let property = next_arg(tokens, pos, "p.")?;
        return parse_custom_property_selector(property, tokens, pos);
    }
    if lower.starts_with("p.") && lower.len() > 2 {
        let property = tok[2..].to_string();
        *pos += 1;
        return parse_custom_property_selector(property, tokens, pos);
    }

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
        "state" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "state")?;
            Ok(Selector::State(parse_state_arg(&arg)?))
        }
        "bonded" => {
            *pos += 1;
            Ok(Selector::Bonded)
        }
        "donors" | "don." | "hbd." => {
            *pos += 1;
            Ok(Selector::Donors)
        }
        "acceptors" | "acc." | "hba." => {
            *pos += 1;
            Ok(Selector::Acceptors)
        }
        "delocalized" | "deloc." => {
            *pos += 1;
            Ok(Selector::Delocalized)
        }
        "flag" | "f;" | "f." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "flag")?;
            Ok(Selector::Flag(parse_flag_number(&arg)?))
        }
        "fixed" | "fxd." => {
            *pos += 1;
            Ok(Selector::Flag(3))
        }
        "restrained" | "rst." => {
            *pos += 1;
            Ok(Selector::Flag(2))
        }
        "masked" | "msk." => {
            *pos += 1;
            Ok(Selector::Masked)
        }
        "protected" => {
            *pos += 1;
            Ok(Selector::Protected)
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
        "segment" | "segid" | "segi" | "s;" | "s." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "segi")?;
            Ok(Selector::Segi(arg))
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
        "pepseq" | "ps." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "pepseq")?;
            Ok(Selector::Pepseq(arg.to_uppercase()))
        }
        "text_type" | "tt;" | "tt." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "text_type")?;
            Ok(Selector::TextType(arg))
        }
        "numeric_type" | "nt;" | "nt." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "numeric_type")?;
            Ok(Selector::NumericType(parse_numeric_range_list(
                &arg,
                "numeric_type",
            )?))
        }
        "custom" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "custom")?;
            Ok(Selector::Custom(arg))
        }
        "label" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "label")?;
            Ok(Selector::Label(arg))
        }
        "stereo" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "stereo")?;
            Ok(Selector::Stereo(arg))
        }
        "serial" | "id" | "ID" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "serial")?;
            parse_serial_range(&arg)
        }
        "index" | "idx." => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "index")?;
            parse_index_range(&arg)
        }
        "rank" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "rank")?;
            parse_rank_range(&arg)
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
        "cartoon_color" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "cartoon_color")?;
            Ok(Selector::CartoonColor(arg))
        }
        "ribbon_color" => {
            *pos += 1;
            let arg = next_arg(tokens, pos, "ribbon_color")?;
            Ok(Selector::RibbonColor(arg))
        }
        "b" => {
            *pos += 1;
            parse_property_selector(AtomProperty::BFactor, tokens, pos, "b")
        }
        "q" => {
            *pos += 1;
            parse_property_selector(AtomProperty::Occupancy, tokens, pos, "q")
        }
        "formal_charge" | "fc;" | "fc." => {
            *pos += 1;
            parse_property_selector(AtomProperty::FormalCharge, tokens, pos, "formal_charge")
        }
        "partial_charge" | "pc;" | "pc." => {
            *pos += 1;
            parse_property_selector(AtomProperty::PartialCharge, tokens, pos, "partial_charge")
        }
        "vdw" => {
            *pos += 1;
            parse_property_selector(AtomProperty::Vdw, tokens, pos, "vdw")
        }
        "elec_radius" => {
            *pos += 1;
            parse_property_selector(AtomProperty::ElecRadius, tokens, pos, "elec_radius")
        }
        "cartoon" => {
            *pos += 1;
            parse_property_selector(AtomProperty::Cartoon, tokens, pos, "cartoon")
        }
        "geom" => {
            *pos += 1;
            parse_property_selector(AtomProperty::Geom, tokens, pos, "geom")
        }
        "valence" => {
            *pos += 1;
            parse_property_selector(AtomProperty::Valence, tokens, pos, "valence")
        }
        "reps" => {
            *pos += 1;
            parse_property_selector(AtomProperty::Reps, tokens, pos, "reps")
        }
        "protons" => {
            *pos += 1;
            parse_property_selector(AtomProperty::Protons, tokens, pos, "protons")
        }
        "flags" => {
            *pos += 1;
            parse_property_selector(AtomProperty::Flags, tokens, pos, "flags")
        }
        "explicit_degree" => {
            *pos += 1;
            parse_property_selector(AtomProperty::ExplicitDegree, tokens, pos, "explicit_degree")
        }
        "explicit_valence" => {
            *pos += 1;
            parse_property_selector(
                AtomProperty::ExplicitValence,
                tokens,
                pos,
                "explicit_valence",
            )
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
        "gap" => {
            *pos += 1;
            let distance = parse_distance_arg(tokens, pos, "gap")?;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Gap(distance, Box::new(inner)))
        }
        "byresidue" | "byresi" | "byres" | "br;" | "br." | "b;" => {
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
        "bysegment" | "byseg" | "bysegi" | "bs." => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Bysegment(Box::new(inner)))
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
        "byfragment" | "byfrag" | "bf." => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Bymolecule(Box::new(inner)))
        }
        "byring" => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Byring(Box::new(inner)))
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
        "neighbor" | "nbr" | "nbr;" | "nbr." => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::Neighbor(Box::new(inner)))
        }
        "bound_to" | "bto." => {
            *pos += 1;
            let inner = parse_unary_inner(tokens, pos)?;
            Ok(Selector::BoundTo(Box::new(inner)))
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
        _ => {
            *pos += 1;
            Ok(Selector::Identifier(tok.to_string()))
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

fn is_and_operator(token: &str) -> bool {
    token == "&" || token.eq_ignore_ascii_case("and")
}

fn is_or_operator(token: &str) -> bool {
    token == "|" || token == "+" || token.eq_ignore_ascii_case("or")
}

fn is_subtract_operator(token: &str) -> bool {
    token == "-"
}

fn is_in_operator(token: &str) -> bool {
    token.eq_ignore_ascii_case("in")
}

fn is_like_operator(token: &str) -> bool {
    token.eq_ignore_ascii_case("like") || token == "l;" || token == "l."
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

fn parse_custom_property_selector(
    property: String,
    tokens: &[String],
    pos: &mut usize,
) -> Result<Selector, String> {
    if property.trim().is_empty() {
        return Err("'p.' requires a property name".to_string());
    }

    let op_token = next_arg(tokens, pos, "p.")?;
    let value_token = next_arg(tokens, pos, "p.")?;
    let op = parse_custom_property_op(&op_token).ok_or_else(|| {
        format!(
            "Invalid comparison operator for 'p.{}': '{}'",
            property, op_token
        )
    })?;

    if !matches!(op, CustomPropertyOp::In) {
        value_token.parse::<f32>().map_err(|_| {
            format!(
                "Invalid numeric comparison value for 'p.{}': '{}'",
                property, value_token
            )
        })?;
    }

    Ok(Selector::CustomProperty(property, op, value_token))
}

fn parse_distance_arg(tokens: &[String], pos: &mut usize, keyword: &str) -> Result<f32, String> {
    let distance_token = next_arg(tokens, pos, keyword)?;
    distance_token
        .parse()
        .map_err(|_| format!("Invalid distance '{}'", distance_token))
}

fn parse_unary_inner(tokens: &[String], pos: &mut usize) -> Result<Selector, String> {
    if *pos < tokens.len() && tokens[*pos] == "of" {
        *pos += 1;
        if *pos >= tokens.len() {
            return Err("Expected selection after 'of'".to_string());
        }
    }

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

fn parse_custom_property_op(s: &str) -> Option<CustomPropertyOp> {
    match s {
        ">" => Some(CustomPropertyOp::Greater),
        "<" => Some(CustomPropertyOp::Less),
        "=" | "==" => Some(CustomPropertyOp::Equal),
        ">=" => Some(CustomPropertyOp::GreaterEqual),
        "<=" => Some(CustomPropertyOp::LessEqual),
        _ if s.eq_ignore_ascii_case("in") => Some(CustomPropertyOp::In),
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
        Ok(Selector::ResiList(parse_resi_range_list(s, "resi")?))
    } else {
        let (start, end, ins_lo, ins_hi) = parse_resi_range_single(s, "resi")?;
        Ok(Selector::Resi(start, end, ins_lo, ins_hi))
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

fn parse_index_range(s: &str) -> Result<Selector, String> {
    if s.contains('+') {
        Ok(Selector::IndexList(parse_numeric_range_list(s, "index")?))
    } else {
        let (start, end) = parse_numeric_range(s, "index")?;
        Ok(Selector::Index(start, end))
    }
}

fn parse_rank_range(s: &str) -> Result<Selector, String> {
    if s.contains('+') {
        Ok(Selector::RankList(parse_numeric_range_list(s, "rank")?))
    } else {
        let (start, end) = parse_numeric_range(s, "rank")?;
        Ok(Selector::Rank(start, end))
    }
}

fn parse_flag_number(s: &str) -> Result<u8, String> {
    let flag: u8 = s.parse().map_err(|_| format!("Invalid flag: '{s}'"))?;
    if flag <= 31 {
        Ok(flag)
    } else {
        Err(format!("Invalid flag: '{s}'"))
    }
}

fn parse_state_arg(s: &str) -> Result<isize, String> {
    let state: isize = s.parse().map_err(|_| format!("Invalid state: '{s}'"))?;
    if state == -1 || state >= 1 {
        Ok(state)
    } else {
        Err(format!("Invalid state: '{s}'"))
    }
}

fn parse_numeric_range_list(s: &str, label: &str) -> Result<Vec<(i32, i32)>, String> {
    s.split('+')
        .map(|part| parse_numeric_range(part.trim(), label))
        .collect()
}

fn parse_resi_range_list(
    s: &str,
    label: &str,
) -> Result<Vec<(i32, i32, Option<char>, Option<char>)>, String> {
    s.split('+')
        .map(|part| parse_resi_range_single(part.trim(), label))
        .collect()
}

fn parse_resi_range_single(
    s: &str,
    label: &str,
) -> Result<(i32, i32, Option<char>, Option<char>), String> {
    let s = s.trim();

    if let Some((start, end)) = split_resi_range(s, "-") {
        let (start_num, start_ins) = parse_resi_end(start, label)?;
        let (end_num, end_ins) = parse_resi_end(end, label)?;
        return Ok((start_num, end_num, start_ins, end_ins));
    }

    if let Some((start, end)) = split_resi_range(s, ":") {
        let (start_num, start_ins) = parse_resi_end(start, label)?;
        let (end_num, end_ins) = parse_resi_end(end, label)?;
        return Ok((start_num, end_num, start_ins, end_ins));
    }

    let (number, ins_code) = parse_resi_end(s, label)?;
    Ok((number, number, ins_code, ins_code))
}

fn parse_resi_end(s: &str, label: &str) -> Result<(i32, Option<char>), String> {
    if s.is_empty() {
        return Err(format!("Invalid {label}: ''"));
    }

    let mut chars = s.chars();
    match chars.next_back() {
        None => Err(format!("Invalid {label}: ''")),
        Some(last) if !last.is_ascii_alphabetic() => parse_resi_number(s, label).map(|n| (n, None)),
        Some(last) => {
            let base = chars.as_str();
            if base.is_empty() {
                Err(format!("Invalid {label}: '{s}'"))
            } else {
                let number = parse_resi_number(base, label)?;
                Ok((number, Some(last.to_ascii_uppercase())))
            }
        }
    }
}

fn split_resi_range<'a>(s: &'a str, delimiter: &str) -> Option<(&'a str, &'a str)> {
    if s.is_empty() {
        return None;
    }

    if s.len() == 1 {
        return None;
    }

    if let Some(idx) = s[1..].find(delimiter) {
        let idx = idx + 1;
        if idx + delimiter.len() >= s.len() {
            return None;
        }
        let before = &s[..idx];
        let after = &s[idx + delimiter.len()..];
        if after.is_empty() {
            return None;
        }
        Some((before, after))
    } else {
        None
    }
}

fn parse_resi_number(s: &str, label: &str) -> Result<i32, String> {
    if s.is_empty() {
        return Err(format!("Invalid {label}: ''"));
    }
    s.parse().map_err(|_| format!("Invalid {label}: '{s}'"))
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
        assert!(matches!(
            parse_selection("chain A:C").unwrap(),
            Selector::ChainPattern(ref s) if s == "A:C"
        ));
    }

    #[test]
    fn test_segi_selection() {
        assert!(matches!(
            parse_selection("segi A1").unwrap(),
            Selector::Segi(ref s) if s == "A1"
        ));
        assert!(matches!(
            parse_selection("segment PRO*").unwrap(),
            Selector::Segi(ref s) if s == "PRO*"
        ));
        assert!(matches!(
            parse_selection("s. A+B").unwrap(),
            Selector::Segi(ref s) if s == "A+B"
        ));
    }

    #[test]
    fn test_resi_range() {
        let sel = parse_selection("resi 10-50").unwrap();
        assert!(matches!(sel, Selector::Resi(10, 50, None, None)));
    }

    #[test]
    fn test_resi_alt_range_and_insertion_code() {
        let sel = parse_selection("resi 2:4").unwrap();
        assert!(matches!(sel, Selector::Resi(2, 4, None, None)));
        assert!(matches!(
            parse_selection("resi 9A").unwrap(),
            Selector::Resi(9, 9, Some('A'), Some('A'))
        ));
        assert!(matches!(
            parse_selection("resi 9A-10A").unwrap(),
            Selector::Resi(9, 10, Some('A'), Some('A'))
        ));
        assert!(matches!(
            parse_selection("resi 9-10A").unwrap(),
            Selector::Resi(9, 10, None, Some('A'))
        ));
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
            parse_selection("pepseq AG").unwrap(),
            Selector::Pepseq(ref s) if s == "AG"
        ));
        assert!(matches!(
            parse_selection("ps. A+G").unwrap(),
            Selector::Pepseq(ref s) if s == "A+G"
        ));
        assert!(matches!(
            parse_selection("text_type CT+HC").unwrap(),
            Selector::TextType(ref s) if s == "CT+HC"
        ));
        assert!(matches!(
            parse_selection("custom ligand*").unwrap(),
            Selector::Custom(ref s) if s == "ligand*"
        ));
        assert!(matches!(
            parse_selection("label active*").unwrap(),
            Selector::Label(ref s) if s == "active*"
        ));
        assert!(matches!(
            parse_selection("stereo R+S").unwrap(),
            Selector::Stereo(ref s) if s == "R+S"
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
                if ranges == &vec![
                    (2, 2, None, None),
                    (12, 12, None, None),
                    (9, 9, None, None),
                    (4, 4, None, None)
                ]
        ));
        assert!(matches!(
            parse_selection("resi 1+3-5").unwrap(),
            Selector::ResiList(ref ranges)
                if ranges == &vec![(1, 1, None, None), (3, 5, None, None)]
        ));
        assert!(matches!(
            parse_selection("index 10+20-22").unwrap(),
            Selector::IndexList(ref ranges) if ranges == &vec![(10, 10), (20, 22)]
        ));
        assert!(matches!(
            parse_selection("numeric_type 10+20-22").unwrap(),
            Selector::NumericType(ref ranges) if ranges == &vec![(10, 10), (20, 22)]
        ));
    }

    #[test]
    fn test_selector_aliases() {
        assert!(matches!(parse_selection("*").unwrap(), Selector::All));
        assert!(matches!(
            parse_selection("residue 10").unwrap(),
            Selector::Resi(10, 10, None, None)
        ));
        assert!(matches!(
            parse_selection("i; 10").unwrap(),
            Selector::Resi(10, 10, None, None)
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
            parse_selection("tt. CT").unwrap(),
            Selector::TextType(ref s) if s == "CT"
        ));
        assert!(matches!(
            parse_selection("nt. 42").unwrap(),
            Selector::NumericType(ref ranges) if ranges == &vec![(42, 42)]
        ));
        assert!(matches!(
            parse_selection("nt; 42").unwrap(),
            Selector::NumericType(ref ranges) if ranges == &vec![(42, 42)]
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
        assert!(matches!(
            parse_selection("%stored").unwrap(),
            Selector::Named(ref name) if name == "stored"
        ));
        assert!(matches!(
            parse_selection("% stored").unwrap(),
            Selector::Named(ref name) if name == "stored"
        ));
        assert!(matches!(
            parse_selection("stored").unwrap(),
            Selector::Identifier(ref name) if name == "stored"
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
            parse_selection("state 2").unwrap(),
            Selector::State(2)
        ));
        assert!(matches!(
            parse_selection("state -1").unwrap(),
            Selector::State(-1)
        ));
        assert!(parse_selection("state 0").is_err());
        assert!(matches!(
            parse_selection("bonded").unwrap(),
            Selector::Bonded
        ));
        assert!(matches!(
            parse_selection("donors").unwrap(),
            Selector::Donors
        ));
        assert!(matches!(parse_selection("don.").unwrap(), Selector::Donors));
        assert!(matches!(parse_selection("hbd.").unwrap(), Selector::Donors));
        assert!(matches!(
            parse_selection("acceptors").unwrap(),
            Selector::Acceptors
        ));
        assert!(matches!(
            parse_selection("acc.").unwrap(),
            Selector::Acceptors
        ));
        assert!(matches!(
            parse_selection("hba.").unwrap(),
            Selector::Acceptors
        ));
        assert!(matches!(
            parse_selection("delocalized").unwrap(),
            Selector::Delocalized
        ));
        assert!(matches!(
            parse_selection("deloc.").unwrap(),
            Selector::Delocalized
        ));
        assert!(matches!(
            parse_selection("flag 25").unwrap(),
            Selector::Flag(25)
        ));
        assert!(matches!(
            parse_selection("f. 31").unwrap(),
            Selector::Flag(31)
        ));
        assert!(matches!(
            parse_selection("f; 0").unwrap(),
            Selector::Flag(0)
        ));
        assert!(matches!(
            parse_selection("fixed").unwrap(),
            Selector::Flag(3)
        ));
        assert!(matches!(
            parse_selection("fxd.").unwrap(),
            Selector::Flag(3)
        ));
        assert!(matches!(
            parse_selection("restrained").unwrap(),
            Selector::Flag(2)
        ));
        assert!(matches!(
            parse_selection("rst.").unwrap(),
            Selector::Flag(2)
        ));
        assert!(parse_selection("flag 32").is_err());
        assert!(matches!(
            parse_selection("masked").unwrap(),
            Selector::Masked
        ));
        assert!(matches!(parse_selection("msk.").unwrap(), Selector::Masked));
        assert!(matches!(
            parse_selection("protected").unwrap(),
            Selector::Protected
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
            parse_selection("id 5").unwrap(),
            Selector::Serial(5, 5)
        ));
        assert!(matches!(
            parse_selection("ID 5").unwrap(),
            Selector::Serial(5, 5)
        ));
        assert!(matches!(
            parse_selection("index 10-20").unwrap(),
            Selector::Index(10, 20)
        ));
        assert!(matches!(
            parse_selection("idx. 10-20").unwrap(),
            Selector::Index(10, 20)
        ));
        assert!(matches!(
            parse_selection("idx. 1+3").unwrap(),
            Selector::IndexList(ref ranges) if ranges == &vec![(1, 1), (3, 3)]
        ));
        assert!(matches!(
            parse_selection("rank 10-20").unwrap(),
            Selector::Rank(10, 20)
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
        assert!(matches!(
            parse_selection("cartoon_color red").unwrap(),
            Selector::CartoonColor(ref s) if s == "red"
        ));
        assert!(matches!(
            parse_selection("ribbon_color blue").unwrap(),
            Selector::RibbonColor(ref s) if s == "blue"
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
            parse_selection("formal_charge = -1").unwrap(),
            Selector::Property(AtomProperty::FormalCharge, CompareOp::Equal, v)
                if (v + 1.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("fc. >= 1").unwrap(),
            Selector::Property(AtomProperty::FormalCharge, CompareOp::GreaterEqual, v)
                if (v - 1.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("partial_charge < -0.25").unwrap(),
            Selector::Property(AtomProperty::PartialCharge, CompareOp::Less, v)
                if (v + 0.25).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("pc; > 0.1").unwrap(),
            Selector::Property(AtomProperty::PartialCharge, CompareOp::Greater, v)
                if (v - 0.1).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("vdw <= 1.7").unwrap(),
            Selector::Property(AtomProperty::Vdw, CompareOp::LessEqual, v)
                if (v - 1.7).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("elec_radius > 1.2").unwrap(),
            Selector::Property(AtomProperty::ElecRadius, CompareOp::Greater, v)
                if (v - 1.2).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("cartoon = 2").unwrap(),
            Selector::Property(AtomProperty::Cartoon, CompareOp::Equal, v)
                if (v - 2.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("geom >= 3").unwrap(),
            Selector::Property(AtomProperty::Geom, CompareOp::GreaterEqual, v)
                if (v - 3.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("valence < 4").unwrap(),
            Selector::Property(AtomProperty::Valence, CompareOp::Less, v)
                if (v - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("reps = 2").unwrap(),
            Selector::Property(AtomProperty::Reps, CompareOp::Equal, v)
                if (v - 2.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("protons == 6").unwrap(),
            Selector::Property(AtomProperty::Protons, CompareOp::Equal, v)
                if (v - 6.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("flags > 0").unwrap(),
            Selector::Property(AtomProperty::Flags, CompareOp::Greater, v)
                if v.abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("explicit_degree >= 2").unwrap(),
            Selector::Property(AtomProperty::ExplicitDegree, CompareOp::GreaterEqual, v)
                if (v - 2.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("explicit_valence < 4").unwrap(),
            Selector::Property(AtomProperty::ExplicitValence, CompareOp::Less, v)
                if (v - 4.0).abs() < f32::EPSILON
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
        assert!(matches!(
            parse_selection("p.score > 0.5").unwrap(),
            Selector::CustomProperty(ref name, CustomPropertyOp::Greater, ref value)
                if name == "score" && value == "0.5"
        ));
        assert!(matches!(
            parse_selection("p. score <= 1.25").unwrap(),
            Selector::CustomProperty(ref name, CustomPropertyOp::LessEqual, ref value)
                if name == "score" && value == "1.25"
        ));
        assert!(matches!(
            parse_selection("p.kind in ligand*").unwrap(),
            Selector::CustomProperty(ref name, CustomPropertyOp::In, ref value)
                if name == "kind" && value == "ligand*"
        ));
    }

    #[test]
    fn test_within_selector() {
        assert!(matches!(
            parse_selection("within 4 chain A").unwrap(),
            Selector::Within(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("within 4 of chain A").unwrap(),
            Selector::Within(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("w. 4 chain A").unwrap(),
            Selector::Within(d, _) if (d - 4.0).abs() < f32::EPSILON
        ));
        assert!(matches!(
            parse_selection("w. 4 of chain A").unwrap(),
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
            parse_selection("a. 4 of chain A").unwrap(),
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
        assert!(matches!(
            parse_selection("gap 2.5 chain A").unwrap(),
            Selector::Gap(d, _) if (d - 2.5).abs() < f32::EPSILON
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
        assert!(matches!(
            parse_selection("nbr; chain A").unwrap(),
            Selector::Neighbor(_)
        ));
        assert!(matches!(
            parse_selection("nbr. chain A").unwrap(),
            Selector::Neighbor(_)
        ));
        assert!(matches!(
            parse_selection("bound_to chain A").unwrap(),
            Selector::BoundTo(_)
        ));
        assert!(matches!(
            parse_selection("bto. chain A").unwrap(),
            Selector::BoundTo(_)
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
            parse_selection("b; chain A").unwrap(),
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
            parse_selection("bysegment name CA").unwrap(),
            Selector::Bysegment(_)
        ));
        assert!(matches!(
            parse_selection("byseg name CA").unwrap(),
            Selector::Bysegment(_)
        ));
        assert!(matches!(
            parse_selection("bysegi name CA").unwrap(),
            Selector::Bysegment(_)
        ));
        assert!(matches!(
            parse_selection("bs. name CA").unwrap(),
            Selector::Bysegment(_)
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
        assert!(matches!(
            parse_selection("byfragment serial 10").unwrap(),
            Selector::Bymolecule(_)
        ));
        assert!(matches!(
            parse_selection("byfrag serial 10").unwrap(),
            Selector::Bymolecule(_)
        ));
        assert!(matches!(
            parse_selection("bf. serial 10").unwrap(),
            Selector::Bymolecule(_)
        ));
        assert!(matches!(
            parse_selection("byring serial 10").unwrap(),
            Selector::Byring(_)
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

        let sel = parse_selection("chain A + chain B").unwrap();
        assert!(matches!(sel, Selector::Or(_, _)));

        let sel = parse_selection("all - chain B").unwrap();
        assert!(matches!(sel, Selector::And(_, _)));

        let sel = parse_selection("name CA in chain A").unwrap();
        assert!(matches!(sel, Selector::In(_, _)));

        let sel = parse_selection("name CA like chain A").unwrap();
        assert!(matches!(sel, Selector::Like(_, _)));

        let sel = parse_selection("name CA l. chain A").unwrap();
        assert!(matches!(sel, Selector::Like(_, _)));
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
