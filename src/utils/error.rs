use thiserror::Error;
use std::io;

#[derive(Error, Debug)]
pub enum MycoNoteError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    
    #[error("Parse error at line {line}: {message}")]
    ParseError {
        line: usize,
        message: String,
    },
    
    #[error("Invalid GFF format: {0}")]
    InvalidFormat(String),
    
    #[error("Feature error: {0}")]
    FeatureError(String),
    
    #[error("Unsupported file format: {0}")]
    UnsupportedFormat(String),
    
    #[error("Plotting error: {0}")]
    PlottingError(String),

    #[error("External tool error: {0}")]
    ExternalTool(String),

    #[error("JSON serialization error: {0}")]
    JsonError(#[from] serde_json::Error),
}

impl From<plotters::drawing::DrawingAreaErrorKind<plotters_bitmap::BitMapBackendError>> for MycoNoteError {
    fn from(err: plotters::drawing::DrawingAreaErrorKind<plotters_bitmap::BitMapBackendError>) -> Self {
        MycoNoteError::PlottingError(format!("{}", err))
    }
}

pub type Result<T> = std::result::Result<T, MycoNoteError>;
