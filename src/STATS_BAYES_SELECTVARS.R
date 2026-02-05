#/***********************************************************************
# * (C) Copyright Jon K Peck, 2024
# ************************************************************************/

# version 0.1.0

# history
# Jan-2026    Initial version



# helpers
gtxt <- function(...) {
    return(gettext(...,domain="STATS_BAYES_SELECTVARS"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_BAYES_SELECTVARS"))
}

loadmsg = "The R %s package is required but could not be loaded."
tryCatch(suppressWarnings(suppressPackageStartupMessages(library(BAS, warn.conflicts=FALSE))), error=function(e){
    stop(gtxtf(loadmsg,"BAS"), call.=FALSE)
}
)

mylist2env = function(alist) {
    env = new.env()
    lnames = names(alist)
    for (i in 1:length(alist)) {
        assign(lnames[[i]],value = alist[[i]], envir=env)
    }
    return(env)
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = mylist2env(lcl) # makes this list into an environment
    
    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.
        
        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 
        
        if (is.null(msg) || dostop) {
            spssdata.CloseDataConnection()
            lcl$display(inproc)  # display messages and end procedure state

            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any
        

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
                procok = TRUE
            }
        } else {
            procok = inproc
            if (!inproc) {
                procok =tryCatch({
                    spsspkg.StartProcedure(lcl$procname, lcl$omsid)
                    procok = TRUE
                },
                error = function(e) {
                    prockok = FALSE
                }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings and Messages","Warnings", isSplit=FALSE) # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row,
                                               gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory,
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}


casecorrect = function(vlist, vardict, warns) {
    # correct the case of variable names
    # vlist is a list of names, possibly including TO and ALL
    # vardict is a variable dictionary
    # unrecognized names are returned as is as the GetDataFromSPSS api will handle them

    dictnames = vardict["varName",]
    names(dictnames) = tolower(dictnames)
    dictnames['all'] = "all"
    dictnames['to'] = "to"
    correctednames = list()
    for (item in vlist) {
        lcitem = tolower(item)
        itemc = dictnames[[lcitem]]
        if (is.null(itemc)) {
            warns$warn(gtxtf("Invalid variable name: %s", item), dostop=TRUE)
        }
        correctednames = append(correctednames, itemc)
    }
    return(correctednames)
}

procname=gtxt("bayesselectvar")
warningsprocname = gtxt("Bayes Select Vars Notes and Warnings")
omsid="STATSBAYESSELECTVARS"



# main worker
dobayesselectvars<-function(depvar, indvars, forcedvars=NULL, plots=TRUE,
        family="linear", offset=NULL, betaprior="cch", modelprior=NULL,
        priorparams=list()
    ) {
    #DEBUG
    # sink(file="c:/temp/bayesselvarsout.log", type="output")
    # f = file("c:/temp/bayesselvars.log", open="w")
    # sink(file=f, type="message")
    
    domain<-"STATS_BAYES_SELECTVARS"
    setuplocalization(domain)
    warns = Warn(procname=warningsprocname,omsid=omsid)
    
    spsspkg.StartProcedure(gtxt("Bayes Select Variables"),"STATS BAYES SELECTVARS")
    # if (!spsspkg.IsUTF8mode()) {
    #     warns$warn(gtxt("This procedure requires SPSS to be in Unicode mode"), dostop=TRUE)
    # }
    weightvar = spssdictionary.GetWeightVariable()
    if (!is.null(weightvar)) {
        warns$warn(gtxt("Case weight are ignored by this procedure."))
    }
    vardict = spssdictionary.GetDictionaryFromSPSS()

    depvar = unlist(casecorrect(list(depvar), vardict, warns))
    indvars = casecorrect(indvars, vardict, warns)
    forcedvars = casecorrect(forcedvars, vardict, warns)
    allind = unique(c(indvars, forcedvars))

    # todo: allow splits if not saving new variables
    nsplitvars = length(spssdata.GetSplitVariableNames())
    # if (nsplitvars > 0 && newvars) {
    #     warns$warn(gtxt("Split files is not supported by this procedure"), dostop=TRUE)
    # }

    if (length(intersect(depvar, allind)) > 0) {
        warns$warn(gtxt("The dependent variable appears in the independent variable list"),
            dostop=TRUE)
    }
    indvarsplus = paste(c(indvars, forcedvars), collapse="+")
    allvars = c(depvar, allind)
    # check for valid names in R
    tryCatch(
        {
        for (ch in allvars) {
            xxx = str2lang(ch)
        }
        },
        error = function(e) {warns$warn(gtxtf("%s is not a valid R variable name.  Please rename it", ch),
            dostop=TRUE)}
    )
    
    f = paste(depvar, indvarsplus, sep="~", collapse="")
    ncaseslist = list()
    
    warns$warn(gtxtf("Results computed by R package BAS, version:%s", packageVersion("BAS")),
        dostop=FALSE)
    # if (length(forcedvars) > 0) {
    #     warns$warn(gtxtf("Forced variables: %s", paste(forcedvars, collapse=" ")))
    # }
    if (nsplitvars > 0) {
        warns$warn(gtxt("Best model definitions may differ across splits."))
    }

    # get data api requires case match
    while(!spssdata.IsLastSplit()){
        tryCatch(
            {
            dta = spssdata.GetSplitDataFromSPSS(allvars, missingValueToNA=TRUE, factorMode="levels",
                keepUserMissing=FALSE)
    
            },
            error=function(e) {warns$warn(paste(gtxt("error fetching data"), e, sep="\n"), dostop=TRUE)}
        )
    
        dta = dta[complete.cases(dta),]
        # if (is.null(offset) && family != "linear") {
        #     dta['x y z'] = rep(0, nrow(dta))
        #     offset = dta['x y z']
        # }
        gc()
        ncases = nrow(dta)
        niv = length(allvars)  # including the intercept
        if (ncases == 0) {
            warns$warn(gtxt("There are no complete cases in the data"), dostop=TRUE)
        }
        # GLM options
        if (family == "logistic") {
            familyx = binomial(link = "logit")
        } else if (family == "poisson") {
            familyx = poisson(link = "log")
        } else if (family == "gamma") {
            familyx = Gamma(link = 'log')
        }
        
        arglist = list(formula=f, data=dta)
        estimatorx = ifelse(family == "linear", "bas.lm", "bas.glm")
        if (length(forcedvars) > 0) {
            fv =  paste("~", paste(forcedvars,  collapse ="+"))
            arglist = append(arglist, list(include.always=fv))
        }
        if (family == "linear" && !is.null(offset)) {
            warns$warn(gtxt("Offset is not valid for linear models"), dostop=TRUE)
        }
        if (!is.null(offset)) {
            arglist = append(arglist, list(offset=offset))
        }
        if (family != "linear") {
            arglist = append(arglist, list(family=familyx))
        }
        ###save(f, dta, family, estimatorx, allvars,  allind, forcedvars, arglist, file="c:/temp/est.rdata")
        tryCatch(
            {
            res = do.call(estimatorx, arglist)
            ###save(res, file="c:/temp/res.rdata")
            ncaseslist = append(ncaseslist, res$n)
            },
            error = function(e) {warns$warn(paste(gtxt("error estimating equation"), e, sep="\n"), dostop=TRUE)},
            warning = function(w) {warns$warn(paste(gtxt("error estimating equation"), w, sep="\n"), dostop=TRUE)}
        )
    
        displaytables(res, depvar, indvars, forcedvars, method, ncaseslist, vardict, warns)
        if (plots) {
            tryCatch(
                {
                plot(res, lwd=2, which=1:4, pch=16, sub.caption=gtxtf("Family: %s", family))
                image(res)
                },
            error = function(e) {warns$warn(e, dostop=FALSE)}
            )
        }
    }
    spssdata.CloseDataConnection()

    # DEBUG
    # sink(file=NULL, type="output")
    # sink(file=NULL, type="message")
    warns$display(inproc=FALSE)
}



displaytables = function(res, depvar, indvars, forcedvars, method, ncases, vardict, warns) {
    # display all the tables
        s = data.frame(summary(res))
        if (length(res$include.always) > 1) {
            tag =  sprintf("Forced: %s",  paste(res$namesx[res$include.always], collapse=" "))
        } else {
            tag = ""
        }
        x1 = s[1:(nrow(s)-5),]
        # replace 0/1 coding with X, - for better readability
        x1[2:ncol(x1)] = sapply(x1[2:ncol(x1)], function(f) {ifelse(f==1, "X", "-")})
        colnames(x1)[1] = "P(B != 0 | Y)"
        spsspivottable.Display(
            x = x1,
            title = gtxt("Top Models (Up to Five)"),
            templateName = "STATSBESTMODELS",
            isSplit=TRUE,
            caption = gtxtf("Dependent variable: %s\nX = variable included, - = variable excluded\n%s",
                depvar, tag)
         )
         
        thenames = c(row.names(s[1:(nrow(s) - 5),]), 
            gtxt("Bayes Factor vs. Largest Marginal Likelihood"), gtxt("Posterior Probabilities"), gtxt("R-Squared"),
            gtxt("Number of Variables"), gtxt("Log Marginal Likelihood (Log BF vs Intercept Only Model"))
        thenames = thenames[(length(thenames)-4):length(thenames)]
        x2 = s[(nrow(s)-4):nrow(s),]
        colnames(x2)[1] = "P(B != 0 | Y)"
         spsspivottable.Display(
             x = x2[2:ncol(x2)],
             title = gtxt("Model Statistics"),
             rowlabels = thenames,
             templateName ="STATBAYESPROBS",
             isSplit = TRUE
         )
         coefs = data.frame(coef(res)[1:3])
         row.names(coefs) = row.names(x1) 
         row.names(coefs)[1] = gtxt("(Intercept)")
         colnames(coefs)<-c(gtxt("Mean"), gtxt("Std. Dev"), gtxt("p(beta != 0"))
         spsspivottable.Display(
             x = coefs,
             title = gtxt("Posterior Coefficients"),
             templateName = "STATSBAYESCOEFS",
             isSplit = TRUE,
             caption = gtxtf("Number of cases: %s", paste(ncases, collapse=" "))
         )
}

# createvariables = function(res, dta, vardict, idvar, 
#     resvar, predvar, warns) {
#     # add output variables to active dataset
#     
#     if (length(c(resvar, predvar)) == 2 && tolower(resvar) == tolower(predvar)) {
#         warns$warn(gtxt("The same name was given for residuals and predicted values"), dostop=TRUE)
#     }
#     activevars = tolower(vardict['varName',])
#     if ((!is.null(resvar) && tolower(resvar) %in% activevars) || 
#          (!is.null(predvar) && tolower(predvar) %in% activevars)) {
#         warns$warn(gtxt("The residual or prediction variable already exists."), dostop=TRUE)
#     }
#     dictlist = list()
#     # copy id variable properties
#     dictlist[[1]] = vardict[, vardict['varName',] == idvar]
#     preds = data.frame(row.names(dta))
#     if (!is.null(resvar)) {
#         dictlist[[2]] = c(resvar, "Residuals", 0, "F8.2", "scale")
#         preds[[2]] = res$residuals
#     }
#     if (!is.null(predvar)) {
#         dictlist[[length(dictlist)+1]] = c(predvar, "Predicted Values", 0, "F8.2", "scale")
#         preds[predvar] = res$fitted.values
#     }
#     spsspkg.EndProcedure()
#     newdict = spssdictionary.CreateSPSSDictionary(dictlist)
#     csvtospss(paste("tospss", runif(1,0.05,1), sep=""), newdict, preds, idvar)
#     spsspkg.StartProcedure(gtxt("Perm"),"STATS PERM")
#     warns$display(inproc=TRUE)
# }


# csvtospss = function(preddataset, dict, preds, idvar) {
#     # save a temporary csv file and read into SPSS
#     # preddataset is the datgaset name for the prediction data
#     # dict is the spss dictionary object for the prediction data
#     # preds is the data
#     
#     # due to locale and encoding issues, we can't use a simple Submit
#     # to do the Submit, so a temporary file with forced
#     # encoding setting and INSERT is used
#     
#     csvfile = tempfile("csvpred", tmpdir=tempdir(), fileext=".csv")
#     write.csv(preds, file=csvfile, row.names=FALSE)
#     spsscmd = sprintf('
# * Encoding: UTF-8.
#         PRESERVE.
#         SET DECIMAL DOT.
#         GET DATA  /TYPE=TXT
#         /FILE="%s"
#         /ENCODING="UTF8"
#         /DELCASE=LINE
#         /DELIMITERS=","
#         /QUALIFIER=""""
#         /ARRANGEMENT=DELIMITED
#         /FIRSTCASE=2
#         /VARIABLES=', csvfile)
#     
#     varspecs = list()
#     for (v in 1:ncol(dict)) {
#         if (!strsplit(dict[["varFormat", v]], "\\d+") %in% c('A', 'F')) {
#             dict[["varFormat", v]] = "F"
#         }
#     }
#     for (v in 1:ncol(dict)) {
#         varspecs = append(varspecs, paste(dict[["varName", v]], dict[["varFormat", v]], sep=" "))
#     }
#     varspecs = paste(varspecs, collapse="\n")
#     activedataset = getactivedsname()
#     cmd = paste(spsscmd, varspecs, ".\n", sprintf("dataset name %s.", preddataset), collapse="\n")
#     syntemp = tempfile("csvsyn", tmpdir=tempdir(), fileext=".sps")
#     writeLines(cmd, con=syntemp, useBytes=TRUE)
#     spsspkg.Submit(sprintf("INSERT FILE='%s' ENCODING='UTF8'", syntemp))
#     spsspkg.Submit("RESTORE.")
#     spsspkg.Submit(sprintf("DATASET ACTIVATE %s.", activedataset))
#     spsspkg.Submit("EXECUTE")
#     spsspkg.Submit(sprintf("UPDATE /FILE=* /FILE=%s BY %s.", preddataset, idvar))
#     spsspkg.Submit(sprintf("DATASET CLOSE %s", preddataset))
#     unlink(csvfile)
#     unlink(syntemp)
# }
# # subpunct approximates characters invalid in SPSS variable names
# subpunct = "[-’‘%&'()*+,/:;<=>?\\^`{|}~’]"
# fixnames = function(names) {
#     # return list of legal, nonduplicative SPSS variable names for the input list
#     
#     # dta is a list/vector of names to correct
#     # this function may not perfectly match SPSS name rules
#     
#     newnames = c()
#     newnameslc = c()
#     for (name in names) {
#         newname = gsub(subpunct, "_", name)   # eliminate disallowed characters
#         newname = gsub("(^[0-9])", "X_\\1", newname)  # fix names starting with digit
#         newname = gsub("^\\.|\\.$", "_", newname)  # fix names starting or ending with "."
#         # }
#         # ensure that there are no duplicate names
#         # preserve case but compare caselessly
#         basename = newname
#         for (i in 1:1000) {
#             newnamelc = tolower(newname)
#             if (!newnamelc %in% newnameslc) {
#                 break
#             } else {
#                 newname = paste(basename, i, sep="_")
#                 newnamelc = tolower(newname)
#             }
#         }
#         newnames = append(newnames, newname)
#         newnameslc = append(newnameslc, newnamelc)
#     }
# 
#     return(newnames)
# }

# getactivedsname = function() {
#     # There is no api for this
# 
#     ds = spssdata.GetOpenedDataSetList()
#     spsspkg.Submit("DATASET NAME X44074_60093_")  # renames active dataset
#     ds2 = spssdata.GetOpenedDataSetList()
#     diff = setdiff(ds, ds2)  # find out which one changed
#     spsspkg.Submit("DATASET ACTIVATE X44074_60093_")  # reactivate the previously active one
#     cmd = sprintf("DATASET NAME %s", diff)   # and give it back its name
#     spsspkg.Submit(cmd)
#     return(diff)
# }

setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    if (!is.null(fpath)) {
        bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
    }
} 


Run<-function(args){

    cmdname = args[[1]]
    args <- args[[2]]

    # variable keywords are typed as varname instead of existingvarlist in
    # order to allow for case correction of names later, since the data fetching apis are
    # case sensitive

    oobj <- spsspkg.Syntax(templ=list(
        spsspkg.Template("DEPVAR", subc="", ktype="varname", var="depvar", islist=FALSE),
        spsspkg.Template("INDVARS", subc="", ktype="varname", var="indvars", islist=TRUE),
        spsspkg.Template("FORCEDVARS", subc="", ktype="varname", var="forcedvars", islist=TRUE),
        spsspkg.Template("FAMILY", subc="", ktype="str", var="family", islist=FALSE,
            vallist=list("linear", "logistic", "poisson", "gamma")),
        spsspkg.Template("OFFSET", subc="", ktype="varname", var="offset"),
        ###spsspkg.Template("BETAPRIOR", subc="", ktype="str", var="betaprior"),
        ###spsspkg.Template("MODELPRIOR", subc="", ktype="str", var="modelprior", islist=FALSE),
        ###spsspkg.Template("SAMPLINGMETHOD", subc="", ktype="str", var="method", islist=FALSE),
        ###spsspkg.Template("priorparams", subc="", ktype="str", var="priorparams", islist=TRUE),
        
        spsspkg.Template("PLOTS", subc="DISPLAY", ktype="bool", var="plots", islist=FALSE)
        ))

    if ("HELP" %in% attr(args,"names"))
        helper(cmdname)
    else {
        res <- spsspkg.processcmd(oobj, args, "dobayesselectvars")
    }
}


helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }

    if (exists("spsspkg.helper")) {
        assign("helper", spsspkg.helper)
    }
}
