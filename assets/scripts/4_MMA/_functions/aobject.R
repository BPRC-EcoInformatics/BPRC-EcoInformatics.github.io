#setClass("kernel",representation("function",kpar="list"))
#setClass("kernelMatrix",representation("matrix"),prototype=structure(.Data=matrix()))

setClassUnion("listI", c("list","numeric","vector","integer","matrix"))
setClassUnion("output", c("matrix","factor","vector","logical","numeric","list","integer","NULL"))
setClassUnion("input", c("matrix","list"))
setClassUnion("kfunction", c("function","character"))
setClassUnion("mpinput", c("matrix","data.frame","missing"))
setClassUnion("lpinput", c("list","missing"))
setClassUnion("kpinput", c("kernelMatrix","missing"))



setClass("vm", representation(alpha = "listI", ## since setClassUnion is not  working
                              type = "character",
                              kernelf = "kfunction",
                              kpar = "list",
                              xmatrix = "input",
                              ymatrix = "output",
                              fitted = "output",  
                              lev = "vector",
                              nclass = "numeric",
                              error = "vector",
                              cross = "vector",
                              n.action= "ANY",
                              terms = "ANY",
                              kcall = "call"), contains= "VIRTUAL")
                                        #Generic Vector Machine object 


if(!isGeneric("type")){
  if (is.function("type"))
    fun <- type
  else fun <- function(object) standardGeneric("type")
  setGeneric("type", fun)
}
setMethod("type", "vm", function(object) object@type)
setGeneric("type<-", function(x, value) standardGeneric("type<-"))
setReplaceMethod("type", "vm", function(x, value) {
  x@type <- value
  x
})


if(!isGeneric("kernelf")){
  if (is.function("kernelf"))
    fun <- kernelf
  else fun <- function(object) standardGeneric("kernelf")
  setGeneric("kernelf", fun)
}
setMethod("kernelf", "vm", function(object) object@kernelf)
setGeneric("kernelf<-", function(x, value) standardGeneric("kernelf<-"))
setReplaceMethod("kernelf", "vm", function(x, value) {
  x@kernelf <- value
  x
})

if(!isGeneric("kpar")){
  if (is.function("kpar"))
    fun <- kpar
  else fun <- function(object) standardGeneric("kpar")
  setGeneric("kpar", fun)
}
setMethod("kpar", "vm", function(object) object@kpar)
setGeneric("kpar<-", function(x, value) standardGeneric("kpar<-"))
setReplaceMethod("kpar", "vm", function(x, value) {
  x@kpar <- value
  x
})

if(!isGeneric("kcall")){
  if (is.function("kcall"))
    fun <- kcall
  else fun <- function(object) standardGeneric("kcall")
  setGeneric("kcall", fun)
}
setMethod("kcall", "vm", function(object) object@kcall)
setGeneric("kcall<-", function(x, value) standardGeneric("kcall<-"))
setReplaceMethod("kcall", "vm", function(x, value) {
  x@kcall <- value
  x
})

setMethod("terms", "vm", function(x, ...) x@terms)
setGeneric("terms<-", function(x, value) standardGeneric("terms<-"))
setReplaceMethod("terms", "vm", function(x, value) {
  x@terms <- value
  x
})



if(!isGeneric("xmatrix")){
  if (is.function("xmatrix"))
    fun <- xmatrix
  else fun <- function(object) standardGeneric("xmatrix")
  setGeneric("xmatrix", fun)
}
setMethod("xmatrix", "vm", function(object) object@xmatrix)
setGeneric("xmatrix<-", function(x, value) standardGeneric("xmatrix<-"))
setReplaceMethod("xmatrix", "vm", function(x, value) {
  x@xmatrix <- value
  x
})

if(!isGeneric("ymatrix")){
  if (is.function("ymatrix"))
    fun <- ymatrix
  else fun <- function(object) standardGeneric("ymatrix")
  setGeneric("ymatrix", fun)
}
setMethod("ymatrix", "vm", function(object) object@ymatrix)
setGeneric("ymatrix<-", function(x, value) standardGeneric("ymatrix<-"))
setReplaceMethod("ymatrix", "vm", function(x, value) {
  x@ymatrix <- value
  x
})

setMethod("fitted", "vm", function(object, ...) object@fitted)
setGeneric("fitted<-", function(x, value) standardGeneric("fitted<-"))
setReplaceMethod("fitted", "vm", function(x, value) {
  x@fitted <- value
  x
})

if(!isGeneric("lev")){
  if (is.function("lev"))
    fun <- lev
  else fun <- function(object) standardGeneric("lev")
  setGeneric("lev", fun)
}
setMethod("lev", "vm", function(object) object@lev)
setGeneric("lev<-", function(x, value) standardGeneric("lev<-"))
setReplaceMethod("lev", "vm", function(x, value) {
  x@lev <- value
  x
})

#if(!isGeneric("nclass")){
#  if (is.function("nclass"))
#    fun <- nclass
#  else 
#}
fun <- function(object) standardGeneric("nclass")
setGeneric("nclass", fun)
setMethod("nclass", "vm", function(object) object@nclass)
setGeneric("nclass<-", function(x, value) standardGeneric("nclass<-"))
setReplaceMethod("nclass", "vm", function(x, value) {
  x@nclass <- value
  x
})

if(!isGeneric("alpha")){
  if (is.function("alpha"))
    fun <- alpha
  else fun <- function(object) standardGeneric("alpha")
  setGeneric("alpha", fun)
}
setMethod("alpha", "vm", function(object) object@alpha)
setGeneric("alpha<-", function(x, value) standardGeneric("alpha<-"))
setReplaceMethod("alpha", "vm", function(x, value) {
  x@alpha <- value
  x
})

if(!isGeneric("error")){
  if (is.function("error"))
    fun <- error
  else fun <- function(object) standardGeneric("error")
  setGeneric("error", fun)
}
setMethod("error", "vm", function(object) object@error)
setGeneric("error<-", function(x, value) standardGeneric("error<-"))
setReplaceMethod("error", "vm", function(x, value) {
  x@error <- value
  x
})

if(!isGeneric("cross")){
  if (is.function("cross"))
    fun <- cross
  else fun <- function(object) standardGeneric("cross")
  setGeneric("cross", fun)
}
setMethod("cross", "vm", function(object) object@cross)
setGeneric("cross<-", function(x, value) standardGeneric("cross<-"))
setReplaceMethod("cross", "vm", function(x, value) {
  x@cross <- value
  x
})

#if(!isGeneric("n.action")){
#  if (is.function("n.action"))
#    fun <- n.action
#  else 
#}
fun <- function(object) standardGeneric("n.action")
setGeneric("n.action", fun)
setMethod("n.action", "vm", function(object) object@n.action)
setGeneric("n.action<-", function(x, value) standardGeneric("n.action<-"))
setReplaceMethod("n.action", "vm", function(x, value) {
  x@n.action <- value
  x
})




setClass("ksvmMinstep", representation(param = "list",
                                scaling = "ANY",
                                coef = "ANY",
                                alphaindex = "ANY",
                                b = "numeric",
				obj = "vector",
                                SVindex = "vector",
                                nSV = "numeric",
                                prior = "list",
                                prob.model = "list"
                                ), contains="vm")

if(!isGeneric("param")){
  if (is.function("param"))
    fun <- param
  else fun <- function(object) standardGeneric("param")
  setGeneric("param", fun)
}
setMethod("param", "ksvmMinstep", function(object) object@param)
setGeneric("param<-", function(x, value) standardGeneric("param<-"))
setReplaceMethod("param", "ksvmMinstep", function(x, value) {
  x@param <- value
  x
})

if(!isGeneric("scaling")){
  if (is.function("scaling"))
    fun <- scaling
  else fun <- function(object) standardGeneric("scaling")
  setGeneric("scaling", fun)
}
setMethod("scaling", "ksvmMinstep", function(object) object@scaling)
setGeneric("scaling<-", function(x, value) standardGeneric("scaling<-"))
setReplaceMethod("scaling", "ksvmMinstep", function(x, value) {
  x@scaling<- value
  x
})

if(!isGeneric("obj")){
  if (is.function("obj"))
     fun <- obj
  else fun <- function(object) standardGeneric("obj")
  setGeneric("obj", fun)
}
setMethod("obj", "ksvmMinstep", function(object) object@obj)
setGeneric("obj<-", function(x, value) standardGeneric("obj<-"))
setReplaceMethod("obj", "ksvmMinstep", function(x, value) {
   x@obj<- value
   x
})


setMethod("coef", "ksvmMinstep", function(object, ...) object@coef)
setGeneric("coef<-", function(x, value) standardGeneric("coef<-"))
setReplaceMethod("coef", "ksvmMinstep", function(x, value) {
  x@coef <- value
  x
})

if(!isGeneric("alphaindex")){
  if (is.function("alphaindex"))
    fun <- alphaindex
  else fun <- function(object) standardGeneric("alphaindex")
  setGeneric("alphaindex", fun)
}
setMethod("alphaindex", "ksvmMinstep", function(object) object@alphaindex)
setGeneric("alphaindex<-", function(x, value) standardGeneric("alphaindex<-"))
setReplaceMethod("alphaindex", "ksvmMinstep", function(x, value) {
  x@alphaindex <- value
  x
})

if(!isGeneric("b")){
  if (is.function("b"))
    fun <- b
  else fun <- function(object) standardGeneric("b")
  setGeneric("b", fun)
}
setMethod("b", "ksvmMinstep", function(object) object@b)
setGeneric("b<-", function(x, value) standardGeneric("b<-"))
setReplaceMethod("b", "ksvmMinstep", function(x, value) {
  x@b <- value
  x
})

if(!isGeneric("SVindex")){
  if (is.function("SVindex"))
    fun <- SVindex
  else fun <- function(object) standardGeneric("SVindex")
  setGeneric("SVindex", fun)
}
setMethod("SVindex", "ksvmMinstep", function(object) object@SVindex)
setGeneric("SVindex<-", function(x, value) standardGeneric("SVindex<-"))
setReplaceMethod("SVindex", "ksvmMinstep", function(x, value) {
  x@SVindex <- value
  x
})

if(!isGeneric("nSV")){
  if (is.function("nSV"))
    fun <- nSV
  else fun <- function(object) standardGeneric("nSV")
  setGeneric("nSV", fun)
}
setMethod("nSV", "ksvmMinstep", function(object) object@nSV)
setGeneric("nSV<-", function(x, value) standardGeneric("nSV<-"))
setReplaceMethod("nSV", "ksvmMinstep", function(x, value) {
  x@nSV <- value
  x
})

if(!isGeneric("prior")){
  if (is.function("prior"))
    fun <- prior
  else fun <- function(object) standardGeneric("prior")
  setGeneric("prior", fun)
}
setMethod("prior", "ksvmMinstep", function(object) object@prior)
setGeneric("prior<-", function(x, value) standardGeneric("prior<-"))
setReplaceMethod("prior", "ksvmMinstep", function(x, value) {
  x@prior <- value
  x
})

if(!isGeneric("prob.model")){
  if (is.function("prob.model"))
    fun <- prob.model
  else fun <- function(object) standardGeneric("prob.model")
  setGeneric("prob.model", fun)
}
setMethod("prob.model", "ksvmMinstep", function(object) object@prob.model)
setGeneric("prob.model<-", function(x, value) standardGeneric("prob.model<-"))
setReplaceMethod("prob.model", "ksvmMinstep", function(x, value) {
  x@prob.model <- value
  x
})
