########################################
### methods only for class 'dimInfo' ###
########################################
#' @aliases get.dimInfo,dimInfo,character-method
#' @rdname get.dimInfo-method
setMethod(f='get.dimInfo', signature=c('dimInfo', 'character'),
	definition=function(object, type) {
		if ( !type %in% c('strInfo', 'dimInfo', 'varName', 'strID', 'posIndex') ) {
			stop("get.dimInfo:: argument 'type' is not valid!\n")
		}
		if ( type == 'strInfo' ) {
			return(object@strInfo)
		}
		if ( type == 'dimInfo' ) {
			return(object@dimInfo)
		}
		if ( type == 'varName' ) {
			return(object@vNames)
		}		
		if ( type == 'strID' ) {
			return(object@strID)
		}
		if ( type == 'posIndex' ) {
			return(object@posIndex)
		}			    			
	}
)

#' @aliases set.dimInfo,dimInfo,character,character-method
#' @rdname set.dimInfo-method
setMethod(f='set.dimInfo', signature=c('dimInfo', 'character', 'character'),
	definition=function(object, type, input) { 
		if ( !type %in% c('strID') ) {
			stop("set.dimInfo:: check argument 'type'!\n")
		}
		
		if ( type == 'strID' ) {
			object@strID <- input	
		}
		validObject(object)
		return(object)
	}
)