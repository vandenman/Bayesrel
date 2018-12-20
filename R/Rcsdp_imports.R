
# # import unexported functions from Rcsdp
# # we import the functions to avoid the print of the control files the original functions does with every call
# get.prob.info <- utils::getFromNamespace("get.prob.info", "Rcsdp")
# prepare.data  <- utils::getFromNamespace("prepare.data", "Rcsdp")
# get.solution  <- utils::getFromNamespace("get.solution", "Rcsdp")
# validate.data <- utils::getFromNamespace("validate.data", "Rcsdp")
write.control.file <- utils::getFromNamespace("write.control.file", "Rcsdp")
