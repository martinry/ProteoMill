

setClass("ora_result", representation(
    gene            = "character",
    database        = "data.table",
    pAdjustMethod   = "character",
    output          = "data.table"
    )
)


setClass("userdata", representation(
    raw            = "data.table",
    main           = "data.table",
    rawidentifiers = "data.table",
    identifiers    = "data.table",
    descriptions   = "data.table",
    deoutput       = "data.table"
    )
)