class DimensionMismatchError(ValueError):
    """Dimension mismatch."""
    pass   

class NotWellDefinedError(ValueError):
    """Needs to be well defined."""
    pass  

class MissingCellError(ValueError):
    """A cell is missed when evaluating top cell."""
    pass  