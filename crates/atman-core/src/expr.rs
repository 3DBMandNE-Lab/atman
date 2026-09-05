//! Covariate expressions for design formulas, anchor lists, and summary
//! columns.
//!
//! Grammar (whitespace ignored):
//!
//! ```text
//! expr    := term (('+' | '-') term)*
//! term    := unary (('*' | '/') unary)*
//! unary   := '-' unary | primary
//! primary := NUMBER | IDENT | '(' expr ')' | FUNC '(' expr ')'
//! FUNC    := 'log10' | 'log2' | 'ln' | 'z'
//! IDENT   := [A-Za-z_][A-Za-z0-9_]*
//! ```
//!
//! `z(...)` standardizes its argument over the rows entering a fit, so it
//! is only meaningful as the outermost node; `Expr::parse` rejects it
//! anywhere else. Evaluation of a missing variable yields `None`, as does any
//! non-finite intermediate.

use crate::contrast::zscore;
use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    Num(f64),
    Var(String),
    Neg(Box<Expr>),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Div(Box<Expr>, Box<Expr>),
    Log10(Box<Expr>),
    Log2(Box<Expr>),
    Ln(Box<Expr>),
    Z(Box<Expr>),
}

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Num(f64),
    Ident(String),
    Plus,
    Minus,
    Star,
    Slash,
    LParen,
    RParen,
}

fn tokenize(text: &str) -> Result<Vec<Token>, String> {
    let chars: Vec<char> = text.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c.is_whitespace() {
            i += 1;
            continue;
        }
        match c {
            '+' => out.push(Token::Plus),
            '-' => out.push(Token::Minus),
            '*' => out.push(Token::Star),
            '/' => out.push(Token::Slash),
            '(' => out.push(Token::LParen),
            ')' => out.push(Token::RParen),
            c if c.is_ascii_digit() || c == '.' => {
                let start = i;
                while i < chars.len() && (chars[i].is_ascii_digit() || chars[i] == '.') {
                    i += 1;
                }
                if i < chars.len() && (chars[i] == 'e' || chars[i] == 'E') {
                    let mut j = i + 1;
                    if j < chars.len() && (chars[j] == '+' || chars[j] == '-') {
                        j += 1;
                    }
                    if j < chars.len() && chars[j].is_ascii_digit() {
                        i = j;
                        while i < chars.len() && chars[i].is_ascii_digit() {
                            i += 1;
                        }
                    }
                }
                let s: String = chars[start..i].iter().collect();
                let v = s
                    .parse::<f64>()
                    .map_err(|_| format!("invalid number {s:?} in expression {text:?}"))?;
                out.push(Token::Num(v));
                continue;
            }
            c if c.is_ascii_alphabetic() || c == '_' => {
                let start = i;
                while i < chars.len() && (chars[i].is_ascii_alphanumeric() || chars[i] == '_') {
                    i += 1;
                }
                out.push(Token::Ident(chars[start..i].iter().collect()));
                continue;
            }
            other => {
                return Err(format!(
                    "unexpected character {other:?} in expression {text:?}"
                ))
            }
        }
        i += 1;
    }
    Ok(out)
}

struct Parser<'a> {
    tokens: Vec<Token>,
    pos: usize,
    text: &'a str,
}

impl Parser<'_> {
    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.pos)
    }
    fn next(&mut self) -> Option<Token> {
        let t = self.tokens.get(self.pos).cloned();
        self.pos += 1;
        t
    }
    fn parse_expr(&mut self) -> Result<Expr, String> {
        let mut lhs = self.parse_term()?;
        loop {
            match self.peek() {
                Some(Token::Plus) => {
                    self.next();
                    let rhs = self.parse_term()?;
                    lhs = Expr::Add(Box::new(lhs), Box::new(rhs));
                }
                Some(Token::Minus) => {
                    self.next();
                    let rhs = self.parse_term()?;
                    lhs = Expr::Sub(Box::new(lhs), Box::new(rhs));
                }
                _ => return Ok(lhs),
            }
        }
    }
    fn parse_term(&mut self) -> Result<Expr, String> {
        let mut lhs = self.parse_unary()?;
        loop {
            match self.peek() {
                Some(Token::Star) => {
                    self.next();
                    let rhs = self.parse_unary()?;
                    lhs = Expr::Mul(Box::new(lhs), Box::new(rhs));
                }
                Some(Token::Slash) => {
                    self.next();
                    let rhs = self.parse_unary()?;
                    lhs = Expr::Div(Box::new(lhs), Box::new(rhs));
                }
                _ => return Ok(lhs),
            }
        }
    }
    fn parse_unary(&mut self) -> Result<Expr, String> {
        if let Some(Token::Minus) = self.peek() {
            self.next();
            let inner = self.parse_unary()?;
            return Ok(Expr::Neg(Box::new(inner)));
        }
        self.parse_primary()
    }
    fn parse_primary(&mut self) -> Result<Expr, String> {
        match self.next() {
            Some(Token::Num(v)) => Ok(Expr::Num(v)),
            Some(Token::LParen) => {
                let e = self.parse_expr()?;
                match self.next() {
                    Some(Token::RParen) => Ok(e),
                    _ => Err(format!("missing ')' in expression {:?}", self.text)),
                }
            }
            Some(Token::Ident(name)) => {
                if let Some(Token::LParen) = self.peek() {
                    self.next();
                    let arg = self.parse_expr()?;
                    match self.next() {
                        Some(Token::RParen) => {}
                        _ => {
                            return Err(format!(
                                "missing ')' after {name}( in expression {:?}",
                                self.text
                            ))
                        }
                    }
                    match name.as_str() {
                        "log10" => Ok(Expr::Log10(Box::new(arg))),
                        "log2" => Ok(Expr::Log2(Box::new(arg))),
                        "ln" => Ok(Expr::Ln(Box::new(arg))),
                        "z" => Ok(Expr::Z(Box::new(arg))),
                        other => Err(format!(
                            "unknown function {other:?} in expression {:?}; supported: log10, log2, ln, z",
                            self.text
                        )),
                    }
                } else {
                    Ok(Expr::Var(name))
                }
            }
            Some(tok) => Err(format!(
                "unexpected token {tok:?} in expression {:?}",
                self.text
            )),
            None => Err(format!("unexpected end of expression {:?}", self.text)),
        }
    }
}

impl Expr {
    /// Parse an expression; `z(...)` is accepted only as the outermost node.
    pub fn parse(text: &str) -> Result<Expr, String> {
        let trimmed = text.trim();
        if trimmed.is_empty() {
            return Err("empty expression".to_string());
        }
        let tokens = tokenize(trimmed)?;
        let mut parser = Parser {
            tokens,
            pos: 0,
            text: trimmed,
        };
        let expr = parser.parse_expr()?;
        if parser.pos != parser.tokens.len() {
            return Err(format!("trailing input in expression {trimmed:?}"));
        }
        if expr.inner().contains_z() {
            return Err(format!(
                "z(...) may only be the outermost function in {trimmed:?}"
            ));
        }
        Ok(expr)
    }

    fn contains_z(&self) -> bool {
        match self {
            Expr::Z(_) => true,
            Expr::Num(_) | Expr::Var(_) => false,
            Expr::Neg(a) | Expr::Log10(a) | Expr::Log2(a) | Expr::Ln(a) => a.contains_z(),
            Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                a.contains_z() || b.contains_z()
            }
        }
    }

    /// `Some(name)` when the expression is exactly one identifier.
    pub fn bare_identifier(&self) -> Option<&str> {
        match self {
            Expr::Var(v) => Some(v.as_str()),
            _ => None,
        }
    }

    /// Sorted, de-duplicated identifiers referenced anywhere in the expression.
    pub fn variables(&self) -> Vec<String> {
        let mut set = BTreeSet::new();
        self.collect_vars(&mut set);
        set.into_iter().collect()
    }

    fn collect_vars(&self, out: &mut BTreeSet<String>) {
        match self {
            Expr::Num(_) => {}
            Expr::Var(v) => {
                out.insert(v.clone());
            }
            Expr::Neg(a) | Expr::Log10(a) | Expr::Log2(a) | Expr::Ln(a) | Expr::Z(a) => {
                a.collect_vars(out)
            }
            Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                a.collect_vars(out);
                b.collect_vars(out);
            }
        }
    }

    /// True when the outermost node is `z(...)`.
    pub fn is_standardized(&self) -> bool {
        matches!(self, Expr::Z(_))
    }

    /// The expression inside an outermost `z(...)`, or `self`.
    pub fn inner(&self) -> &Expr {
        match self {
            Expr::Z(a) => a,
            other => other,
        }
    }

    /// Row-level evaluation. A `Z` node evaluates its argument (the caller
    /// standardizes across rows via [`evaluate_column`]).
    pub fn eval(&self, lookup: &dyn Fn(&str) -> Option<f64>) -> Option<f64> {
        let v = match self {
            Expr::Num(v) => *v,
            Expr::Var(name) => lookup(name)?,
            Expr::Neg(a) => -a.eval(lookup)?,
            Expr::Add(a, b) => a.eval(lookup)? + b.eval(lookup)?,
            Expr::Sub(a, b) => a.eval(lookup)? - b.eval(lookup)?,
            Expr::Mul(a, b) => a.eval(lookup)? * b.eval(lookup)?,
            Expr::Div(a, b) => a.eval(lookup)? / b.eval(lookup)?,
            Expr::Log10(a) => a.eval(lookup)?.log10(),
            Expr::Log2(a) => a.eval(lookup)?.log2(),
            Expr::Ln(a) => a.eval(lookup)?.ln(),
            Expr::Z(a) => a.eval(lookup)?,
        };
        if v.is_finite() {
            Some(v)
        } else {
            None
        }
    }
}

/// Evaluate `expr` for rows `0..n`, standardizing over the finite rows when
/// the expression is `z(...)`.
pub fn evaluate_column(
    expr: &Expr,
    n: usize,
    lookup: &dyn Fn(usize, &str) -> Option<f64>,
) -> Vec<Option<f64>> {
    let raw: Vec<Option<f64>> = (0..n)
        .map(|i| expr.inner().eval(&|name| lookup(i, name)))
        .collect();
    if expr.is_standardized() {
        zscore(&raw)
    } else {
        raw
    }
}

/// Split `text` on `sep` occurrences that are outside parentheses; parts are
/// trimmed and empty parts dropped.
pub fn split_top_level(text: &str, sep: char) -> Vec<String> {
    let mut out = Vec::new();
    let mut depth = 0i32;
    let mut current = String::new();
    for c in text.chars() {
        match c {
            '(' => {
                depth += 1;
                current.push(c);
            }
            ')' => {
                depth -= 1;
                current.push(c);
            }
            c if c == sep && depth == 0 => {
                let part = current.trim().to_string();
                if !part.is_empty() {
                    out.push(part);
                }
                current.clear();
            }
            _ => current.push(c),
        }
    }
    let part = current.trim().to_string();
    if !part.is_empty() {
        out.push(part);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lookup(name: &str) -> Option<f64> {
        match name {
            "QAlb" => Some(100.0),
            "leukocyte_count" => Some(9.0),
            "age" => Some(50.0),
            _ => None,
        }
    }

    #[test]
    fn parses_bare_identifier() {
        let e = Expr::parse("age").unwrap();
        assert_eq!(e.bare_identifier(), Some("age"));
        assert_eq!(e.variables(), vec!["age".to_string()]);
        assert_eq!(e.eval(&lookup), Some(50.0));
    }

    #[test]
    fn evaluates_functions_and_arithmetic() {
        assert_eq!(Expr::parse("log10(QAlb)").unwrap().eval(&lookup), Some(2.0));
        assert!(
            (Expr::parse("log10(leukocyte_count + 1)")
                .unwrap()
                .eval(&lookup)
                .unwrap()
                - 1.0)
                .abs()
                < 1e-12
        );
        assert!(
            (Expr::parse("log2(QAlb / 25)")
                .unwrap()
                .eval(&lookup)
                .unwrap()
                - 2.0)
                .abs()
                < 1e-12
        );
        assert!(
            (Expr::parse("ln(QAlb) * 2 - age")
                .unwrap()
                .eval(&lookup)
                .unwrap()
                - (2.0 * 100f64.ln() - 50.0))
                .abs()
                < 1e-12
        );
        assert_eq!(Expr::parse("-age + 60").unwrap().eval(&lookup), Some(10.0));
        assert_eq!(Expr::parse("QAlb / missing").unwrap().eval(&lookup), None);
        assert_eq!(Expr::parse("2.5e1 + 1").unwrap().eval(&lookup), Some(26.0));
    }

    #[test]
    fn z_is_only_allowed_outermost() {
        let e = Expr::parse("z(age)").unwrap();
        assert!(e.is_standardized());
        assert_eq!(e.inner().bare_identifier(), Some("age"));
        assert!(Expr::parse("log10(z(age))").is_err());
        assert!(Expr::parse("z(age) + 1").is_err());
    }

    #[test]
    fn rejects_malformed_input() {
        assert!(Expr::parse("").is_err());
        assert!(Expr::parse("log10(").is_err());
        assert!(Expr::parse("age +").is_err());
        assert!(Expr::parse("sqrt(age)").is_err());
        assert!(Expr::parse("age age").is_err());
    }

    #[test]
    fn variables_are_sorted_and_unique() {
        let e = Expr::parse("log10(QIgG / QAlb) + QAlb").unwrap();
        assert_eq!(e.variables(), vec!["QAlb".to_string(), "QIgG".to_string()]);
    }

    #[test]
    fn evaluate_column_standardizes_over_finite_rows() {
        let ages = [Some(10.0), Some(20.0), None, Some(30.0)];
        let e = Expr::parse("z(age)").unwrap();
        let col = evaluate_column(&e, 4, &|i, name| if name == "age" { ages[i] } else { None });
        assert!((col[0].unwrap() + 1.0).abs() < 1e-12);
        assert!((col[1].unwrap()).abs() < 1e-12);
        assert!(col[2].is_none());
        assert!((col[3].unwrap() - 1.0).abs() < 1e-12);
        let raw = evaluate_column(&Expr::parse("age").unwrap(), 4, &|i, name| {
            if name == "age" {
                ages[i]
            } else {
                None
            }
        });
        assert_eq!(raw, vec![Some(10.0), Some(20.0), None, Some(30.0)]);
    }

    #[test]
    fn split_top_level_respects_parentheses() {
        assert_eq!(
            split_top_level("case + z(age) + log10(leukocyte_count + 1)", '+'),
            vec![
                "case".to_string(),
                "z(age)".to_string(),
                "log10(leukocyte_count + 1)".to_string()
            ]
        );
        assert_eq!(
            split_top_level("log10(QAlb),log10(QIgG/QAlb), age", ','),
            vec!["log10(QAlb)", "log10(QIgG/QAlb)", "age"]
        );
        assert!(split_top_level("", '+').is_empty());
    }
}
