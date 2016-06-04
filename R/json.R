# The following code was lifted from Alex Couture-Beil <rjson_pkg at mofo.ca>'s
# rjson package (version 0.2.15) The C-code was removed and only native R code retained
toJSON <- function(x)
{
    #convert factors to characters
    if( is.factor( x ) == TRUE ) {
        tmp_names <- names( x )
        x = as.character( x )
        names( x ) <- tmp_names
    }

    if( !is.vector(x) && !is.null(x) && !is.list(x) ) {
        x <- as.list( x )
        warning("JSON only supports vectors and lists - But I'll try anyways")
    }
    
    if( is.null(x) )
        return( "null" )

    #treat named vectors as lists
    if( is.null( names( x ) ) == FALSE ) {
        x <- as.list( x )
    }
    
    #named lists only
    if( is.list(x) && !is.null(names(x)) ) {
        if( any(duplicated(names(x))) )
            stop( "A JSON list must have unique names" );
        str = "{"
        first_elem = TRUE
        for( n in names(x) ) {
            if( first_elem )
                first_elem = FALSE
            else
                str = paste(str, ',', sep="")
            str = paste(str, deparse(n), ":", toJSON(x[[n]]), sep="")
        }
        str = paste( str, "}", sep="" )
        return( str )
    }
    
    #treat lists without names as JSON array
    if( length(x) != 1 || is.list(x) ) {
        if( !is.null(names(x)) )
            return( toJSON(as.list(x)) ) #vector with names - treat as JSON list
        str = "["
        first_elem = TRUE
        for( val in x ) {
            if( first_elem )
                first_elem = FALSE
            else
                str = paste(str, ',', sep="")
            str = paste(str, toJSON(val), sep="")
        }
        str = paste( str, "]", sep="" )
        return( str )
    }

    if( is.nan(x) )
        return( "\"NaN\"" )

    if( is.na(x) )
        return( "\"NA\"" )

    if( is.infinite(x) )
        return( ifelse( x == Inf, "\"Inf\"", "\"-Inf\"" ) )
    
    if( is.logical(x) )
        return( ifelse(x, "true", "false") )
    
    if( is.character(x) )
        return( gsub("\\/", "\\\\/", deparse(x)) )
    
    if( is.numeric(x) )
        return( as.character(x) )
    
    stop( "shouldnt make it here - unhandled type not caught" )
}

#create an object, which can be used to parse JSON data spanning multiple buffers
#it will be able to pull out multiple objects.. e.g: "[5][2,1]" is two different
# JSON objects - it can be called twice to get both items
newJSONParser <- function()
{
    buffer <- c()
    return(	list(
        "addData" = function( buf ) { 
            chars = strsplit(buf, "")[[1]]
            for( ch in chars )
                buffer[ length(buffer) + 1 ]  <<- ch
        },
        "getObject" = function()
        {
            tmp <- .parseValue( buffer, 1)
            if( is.null( tmp$incomplete ) == FALSE )
                return( NULL )

            if( tmp$size > length(buffer) )
                buffer <<- c()
            else
                buffer <<- buffer[ tmp$size : length(buffer) ]

            return( tmp$val )
        }
    ) )
}


fromJSON <- function( json_str, file, unexpected.escape = "error" )
{
    if( missing( json_str ) ) {
        if( missing( file ) )
            stop( "either json_str or file must be supplied to fromJSON")
        json_str <- paste(readLines( file, warn=FALSE ),collapse="")
    } else {
        if( missing( file ) == FALSE ) {
            stop( "only one of json_str or file must be supplied to fromJSON")
        }
    }
    return( .fromJSON_R( json_str ) )
}

.fromJSON_R <- function( json_str )
{
    if( !is.character(json_str) )
        stop( "JSON objects must be a character string" )
    chars = strsplit(json_str, "")[[1]]
    tmp <- .parseValue( chars, 1)
    if( is.null( tmp$incomplete ) )
        return( tmp$val )
    else
        return( NULL )
}

.parseValue <- function( chars, i )
{
    if( i > length( chars ) )
        return( list( "incomplete" = TRUE ) )
    
    #ignore whitespace
    while( chars[i] == " " || chars[i] == "\t" || chars[i] == "\n" ) {
        i = i + 1
        if( i > length( chars ) )
            return( list( "incomplete" = TRUE ) )
    }

    ch = chars[i]
    if( ch == "{" ) {
        return( .parseObj( chars, i ) )
    }
    if( ch == "[" ) {
        return( .parseArray( chars, i ) )
    }
    if( ch == "\"" ) {
        return( .parseString( chars, i ) )
    }
    if( any(grep("[0-9\\-]", ch)) ) {
        return( .parseNumber( chars, i ) )
    }
    if( ch == "t" ) {
        return( .parseTrue( chars, i ) )
    }
    if( ch == "f" ) {
        return( .parseFalse( chars, i ) )
    }
    if( ch == "n" ) {
        return( .parseNull( chars, i ) )
    }
    #stop("shouldnt reach end of parseValue")
    
    err <- paste( "unexpected data:", paste( chars[ i:length(chars)], collapse = "" ) )
    stop( err )
}

.parseObj <- function( chars, i )
{
    obj <- list()
    if( chars[i] != "{" ) stop("error - no openning tag")
    i = i + 1
    if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
    
    first_pass <- TRUE
    while( TRUE ) {
	
        #ignore whitespace
        while( chars[i] == " " || chars[i] == "\t" || chars[i] == "\n" ) {
            i = i + 1
            if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        }

        
        #look out for empty lists
        if( chars[i] == "}" && first_pass == TRUE ) {
            i = i + 1
            break
        }
        first_pass <- FALSE
        
        #get key
        str = .parseString( chars, i )
        if( is.null( str$incomplete ) == FALSE ) return( str )
        key = str$val
        i = str$size
        if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        
        #ignore whitespace
        while( chars[i] == " " || chars[i] == "\t" || chars[i] == "\n" ) {
            i = i + 1
            if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        }
        
        #verify seperater
        if( chars[i] != ":" ) stop("error - no seperator")
        i = i + 1
        if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )

        
        #ignore whitespace
        while( chars[i] == " " || chars[i] == "\t" || chars[i] == "\n" ) {
            i = i + 1
            if( i > length( chars ) )
                return( list( "incomplete" = TRUE ) )
        }
        
        #get value
        val = .parseValue( chars, i )
        if( is.null( val$incomplete ) == FALSE ) return( val )
        obj[key] <- list(val$val)
        i = val$size
        
        if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
	
        #ignore whitespace
        while( chars[i] == " " || chars[i] == "\t" || chars[i] == "\n" ) {
            i = i + 1
            if( i > length( chars ) )
                return( list( "incomplete" = TRUE ) )
        }
        
        if( chars[i] == "}" ) {
            i = i + 1
            break
        }
        if( chars[i] != "," ) stop("error - no closing tag")
        i = i + 1
        if( i > length( chars ) )
            return( list( "incomplete" = TRUE ) )
    }
    return( list(val=obj, size=i) )
}

.parseArray <- function( chars, i )
{
    useVect <- TRUE
    arr <- list()
    if( chars[i] != "[" ) stop("error - no openning tag")

    i = i + 1
    if( i > length( chars ) )
        return( list( "incomplete" = TRUE ) )

    while( TRUE ) {
        
        #ignore whitespace
        while( chars[i] == " " || chars[i] == "\t" || chars[i] == "\n" ) {
            i = i + 1
            if( i > length( chars ) )
                return( list( "incomplete" = TRUE ) )
        }
	
        #look out for empty arrays
        if( chars[i] == "]" ) { 
            i = i + 1
            useVect <- FALSE #force an empty list instead of NULL (i.e. value = vector("list",0))
            break
        }
        
        #get value
        val = .parseValue( chars, i )
        if( is.null( val$incomplete ) == FALSE ) return( val )
        arr[length(arr)+1] <- list(val$val)
        if( is.list(val$val) || length(val$val) > 1 || is.null(val$val) )
            useVect <- FALSE
        
        i = val$size
        if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        
        #ignore whitespace
        while( chars[i] == " " || chars[i] == "\t" || chars[i] == "\n" ) {
            i = i + 1
            if( i > length( chars ) )
                return( list( "incomplete" = TRUE ) )
        }
        
        if( chars[i] == "]" ) { 
            i = i + 1
            break
        }
        if( chars[i] != "," ) stop("error - no closing tag")
        i = i + 1
        if( i > length( chars ) )
            return( list( "incomplete" = TRUE ) )
    }
    if( useVect )
    	arr <- unlist(arr)
    return( list(val=arr, size=i) )
}

.parseString <- function( chars, i )
{
    str_start = i
    if( chars[i] != "\"") stop("error")
    i = i + 1
    if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
    
    while( TRUE ) {
        while( chars[i] != "\\" && chars[i] != "\"" ) {
            i = i + 1
            if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        }
        if( chars[i] == "\\" ) {
            i = i + 2 #skip the next char
            if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        }
        else
            break
    }
    str_end = i
    i = i + 1
    return(list(
        val=eval(parse( text=paste(chars[str_start:str_end], collapse="") )), 
        size=i ))
}

.parseNumber <- function( chars, i )
{
    str_start = i

    if( chars[i] == "-" )
        i = i + 1
    
    if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
    
    if( chars[i] == "0" ) {
        i = i + 1
        if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        if( any(grep("[1-9]", chars[i])) ) stop("JSON specs don't allow a number like \"012\"")
    } else if( any(grep("[1-9]", chars[i])) ) {
        i = i + 1
        if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        while( any(grep("[0-9]", chars[i])) ) {
            i = i + 1
            if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        }
    } else {
        stop( "doesn't look like a valid JSON number" )
    }
    
    if( chars[i] == "." ) {
        i = i + 1
        if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        while( any(grep("[0-9]", chars[i])) ) {
            i = i + 1
            if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        }
    }
    
    if( chars[i] == "e" || chars[i] == "E" ) {
        i = i + 1
        if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        if( chars[i] == "-" || chars[i] == "+" )
            i = i + 1
        if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        while( any(grep("[0-9]", chars[i])) ) {
            i = i + 1
            if( i > length( chars ) ) return( list( "incomplete" = TRUE ) )
        }
    }
    str_end = i-1
    
    return(list(
        val=eval(parse( text=paste(chars[str_start:str_end], collapse="") )), 
        size=i ))
}

.parseTrue <- function( chars, i )
{
    if( paste(chars[i:(i+3)], collapse="") == "true" )
        return( list(val=TRUE,size=i+4) )
    stop("error parsing true value (maybe the word starts with t but isnt true)")
}

.parseFalse <- function( chars, i )
{
    if( paste(chars[i:(i+4)], collapse="") == "false" )
        return( list(val=FALSE,size=i+5) )
    stop("error parsing false value (maybe the word starts with f but isnt false)")
}

.parseNull <- function( chars, i )
{
    if( paste(chars[i:(i+3)], collapse="") == "null" )
        return( list(val=NULL,size=i+4) )
    stop("error parsing null value (maybe the word starts with n but isnt null)")
}
