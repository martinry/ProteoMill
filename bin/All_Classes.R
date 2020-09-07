

setClass("ora_result", representation(
    gene            = "character",
    database        = "data.table",
    pAdjustMethod   = "character",
    output          = "data.table"
    )
)
