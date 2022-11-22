# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# utilities functions

phaseWithBeagle <- function(geno,
                            out,
                            beaglePath,
                            nthreads,
                            ref = NULL,
                            skipIfExists = F,
                            pedigree = NULL){

  stopifnot(file.exists(beaglePath))
  cmd <- paste0('java -Xmx6g -jar ', beaglePath)

  cmd <- paste0(cmd, ' gt=', geno, ' nthreads=', nthreads)
  if (!is.null(ref)) {
    cmd <- paste0(cmd, ' ref=', ref)
  }

  if (!is.null(pedigree)) {
    cmd <- paste0(cmd, ' ped=', pedigree)
  }

  cmd <- paste0(cmd, ' out=', out)

  outFile <- paste0(tools::file_path_sans_ext(tools::file_path_sans_ext(out)),
                    '.vcf.gz')
  if (!skipIfExists | !file.exists(outFile)) {
    system(cmd)
  }
  outFile
}


writeVCF <- function(pop, file, keepSNP = NULL){

      if (file.exists(file)) {
        stop("`file` should not already exists.")
      }
      ext1 <- tools::file_ext(file)
      ext2 <- tools::file_ext(tools::file_path_sans_ext(file))
      if (ext1 != "gz" || ext2 != "vcf") {
        message('your file extention is not ".vcf.gz"')
      }

      # Fixed region
      data <- pop$inds[[1]]$haplo$SNPinfo$SNPcoord[, c("chr", "physPos", "SNPid")]
      colnames(data) <- c("#CHROM", "POS", "ID")
      data <- data[order(data$POS), ]
      data <- data[order(data$`#CHROM`), ]


      data$REF <- "A" #sample(c("A","C","T","G"), size = nrow(data), replace = TRUE)
      data$ALT <- "C"
      data$QUAL <- "."
      data$FILTER <- "PASS"
      data$INFO <- "."

      if (is.null(keepSNP)) {
        keepSNP <- data$ID
      } else {
        data <- data[data$ID %in% keepSNP,]
      }

      # Genotype region
      data$FORMAT <- "GT"
      gt <- vapply(pop$inds, function(ind) {
        hap <- do.call(cbind, ind$haplo$values)
        hap <- hap[, keepSNP] # filter SNP
        x <- paste(hap[1, ], hap[2, ], sep = "/")
        names(x) <- colnames(hap)
        x
      }, vector(
        mode = "character",
        length = length(keepSNP)
      ))
      gt <- as.data.frame(gt[data$ID, ])
      data <- cbind(data, gt)


      # Meta region
      meta <- paste("##fileformat=VCFv4.3",
        "##source=\"breedSimulatR\", data in this file are simulated.",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        sep = "\n"
      )

      # write file
      f <- gzfile(file, "w")
      writeLines(text = meta, con = f)
      close(f)
      data.table::fwrite(
        x = data,
        file = file,
        append = TRUE,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
      )
}



#' Try to send a notification to the user.
#'
#' It will first try to send a message using `slackr` using the default user's
#' configuration.
#' If it doesn't success, it will try to send an email to the address registered
#' in the `.env` file: "BAYESOPT_MAIL". (postfix should be installed)
#'
#' @param body content of the notification
#'
#' @examples
sendNotification <- function(body = '') {
  tryCatch({
    time <- Sys.time()

    notifBody <- paste0(
      "Notification from `R` at ", time, " on ", Sys.info()["nodename"], ":\n>",
      gsub("\n", "  \n>", body, fixed = T),
      "  \n  \n"
    )

    slackTry <- sendSlack(notifBody)

    if (slackTry != 0) {
      notifBody <- paste0(
        "Notification from `R` at ", time, ":\n\n",
        body,
        "  \n  \n"
      )
      notifBody <- gsub("\n", "\\\n", notifBody, fixed = T)

      email <- git2r::config()$global$user.email

      sendMail(Sys.getenv("BAYESOPT_MAIL"),
               subject = "R notification",
               body = notifBody)
    }

  }, error = function(err) {
    # Error part, evaluated in case of error:
    message("Error detected. Here's the original error message:")
    message(err)
    # Choose a return value in case of error
    return(1)

  }, warning = function(warn) {
    # Warning part, evaluated in case of warning:
    message("Warning detected. Here's the original warning message:")
    message(warn)
    # Choose a return value in case of warning
    return(2)

  })

}




#' Send an email using the linux command "mail" and postfix
#'
#'
#' @param to destination
#' @param subject
#' @param body
#' @param from sender info
#' @param file optional attached files
#'
#' @return NULL
sendMail <- function(to,
                     subject,
                     body,
                     from = paste0(Sys.info()["user"], "-R@", Sys.info()["nodename"]),
                     file = NULL) {

    # check if system can send mails
    x <- system("hash mail", intern = TRUE)
    if ( ! is.null(attributes(x))) {
        stop("Command 'mail' not found. Please install 'postfix' and 'mailutils'.")
    }

    # attached files
    if (! is.null(file)) {
        stopifnot("file not found" = all(file.exists(file)))
        if (sum(file.size(file)) > 10^7 ) {
            warning("Attached files bigger than 10 Mo ! May not work...")
        }
        fileCmd <- paste0(" --content-filename=\"", basename(file),"\" --attach=\"", file, "\"", collapse = " ")
    } else {
        fileCmd <- ""
    }


    tempfile <- tempfile()
    writeLines(body, tempfile)

    command <- paste0("cat ", tempfile, " | mail ",
                      "-s \"", subject, "\" ",
                      "-aFrom:", from,
                      " ", paste(to, collapse = " "),
                      fileCmd)
    # cat(command)
    system(command)
    file.remove(tempfile)
    invisible(NULL)
}


#' Send message to slack
#'
#' @param body body of the message
#' @param channel destination channel of the message
#' @param file attached file
#'
sendSlack <- function(body,
                      channel = NULL,
                      file = NULL) {

  tryCatch({

    # setup slack connexion (need a "~/.slackr" file)
    slackr::slackr_setup()

    if (is.null(channel)) {
      channel <- Sys.getenv("SLACK_CHANNEL")
    }

    slackr::slackr_msg(body, channel = channel)

    if (!is.null(file)) {
      slackr::slackr_upload(filename = file,
                            channels = channel)
    }
    return(0)

  }, error = function(err) {
    # Error part, evaluated in case of error:
    message("Error detected. Here's the original error message:")
    message(err)
    # Choose a return value in case of error
    return(1)

  }, warning = function(warn) {
    # Warning part, evaluated in case of warning:
    message("Warning detected. Here's the original warning message:")
    message(warn)
    # Choose a return value in case of warning
    return(2)

  })


}
