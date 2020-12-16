### MAIN_SCRIPT.R

start_time <- Sys.time()

#################################################################################################
### SET USER PARAMETERS

MAIN_PATH <- "/Users/abrahammathew/Desktop/lyft_automation_project"

RAW_DATA_PATH <- file.path(MAIN_PATH, "raw_data")
CODE_PATH <- file.path(MAIN_PATH, "R")
FINAL_OUTPUT_PATH <- file.path(MAIN_PATH, "final_output")
LOG_OUTPUT_PATH <- file.path(MAIN_PATH, "logs")
SUMMARY_OUTPUT_PATH <- file.path(MAIN_PATH, "summary_output")

USE_THESE_PACKAGES <- c("data.table", "lubridate", "optparse")

FIRST_NAME_COLUMN_NAME <- "first_name"
LAST_NAME_COLUMN_NAME <- "last_name"

SAVE_THESE_COLUMN_NAMES <- c(
        "pickup_date_(local)", "pickup_time_(local)", "first_name", "last_name", 
        "request_address", "destination_address", "distance_(miles)", "duration_(minutes)", 
        "transaction_amount", "passenger_number", "requester_name"
)

TARGET_VARIABLE <- "transaction_amount"
  
TIMESTAMP_SAVE <- format(Sys.time(), "%Y%m%d")

DO_LOGGING <- FALSE

SAVE_METHOD = c("excel")  # "kable", "pdf", "excel", or "csv"

GENERATE_SUMMARY_OUTPUT <- TRUE 

PASSENGER_NUMBER_COLUMN <- "passenger_number"
PICKUP_DATE_COLUMN <- "pickup_date_(local)"
PICKUP_CITY_COLUMN <- "pickup_city"
TRANSACTION_COLUMN <- "transaction_amount"
DISTANCE_COLUMN <- "distance_(miles)"
DURATION_COLUMN <- "duration_(minutes)"
REQUESTER_COLUMN <- "requester_name"

#################################################################################################
### CHECK IF FILE STRUCTURE IS PRESENT 

if(DO_LOGGING){ 
    log_info("Check if File Structure is Present....")
} else {
    message("Check if File Structure is Present....")
}

if(!all(c("raw_data", "R", "final_output", 
         "logs", "summary_output") %in% list.files(MAIN_PATH))){
  
    dir.create(file.path(MAIN_PATH, "raw_data"), showWarnings = FALSE)
    dir.create(file.path(MAIN_PATH, "R"), showWarnings = FALSE)
    dir.create(file.path(MAIN_PATH, "final_output"), showWarnings = FALSE)
    dir.create(file.path(MAIN_PATH, "logs"), showWarnings = FALSE)
    dir.create(file.path(MAIN_PATH, "summary_output"), showWarnings = FALSE)
    
#  if(!all(c("pdf", "excel", "kable", "csv") %in% list.files(FINAL_OUTPUT_PATH))){
#          
#      dir.create(file.path(FINAL_OUTPUT_PATH, "csv"), showWarnings = FALSE)
#      dir.create(file.path(FINAL_OUTPUT_PATH, "excel"), showWarnings = FALSE)
#      dir.create(file.path(FINAL_OUTPUT_PATH, "kable"), showWarnings = FALSE)
#      dir.create(file.path(FINAL_OUTPUT_PATH, "pdf"), showWarnings = FALSE)
#      
#  }
}

#################################################################################################
### PRELIMINARIES

setwd(MAIN_PATH)

options(scipen = 999)
options(tinytex.verbose = TRUE)

library(logger)

log_threshold(TRACE)

log_appender(appender_file(paste("logs/log_file_", TIMESTAMP_SAVE, ".log", sep = "")))

if(DO_LOGGING){ 
    log_info("Script initialized....")
} else {
    message("Script initialized....")
}
  
#################################################################################################
### IMPORT PACKAGES

if(DO_LOGGING){ 
    log_info("Importing packages....")
} else {
    message("Importing packages....")
}
 
new_packages <- USE_THESE_PACKAGES[!(USE_THESE_PACKAGES %in% installed.packages()[,"Package"])]

if(length(new_packages)) {
  install.packages(new_packages)
  sapply(USE_THESE_PACKAGES, require, character.only = TRUE)
} else {
  sapply(USE_THESE_PACKAGES, require, character.only = TRUE)
}

#################################################################################################
### SET FINAL PARAMETERS

if(DO_LOGGING){ 
    log_info("Create List of Parameters....")
} else {
    message("Create List of Parameters....")
}

option_list <- list(

    optparse::make_option(opt_str = c("--MAIN_PATH"), 
                        type = "character", 
                        default = MAIN_PATH, 
                        metavar = "character", 
                        help    = "Root Path for Project"),
    
    optparse::make_option(opt_str = c("--RAW_DATA_PATH"), 
                          type = "character", 
                          default = RAW_DATA_PATH, 
                          metavar = "character", 
                          help    = "Path to Raw Data for Project"),
  
    optparse::make_option(opt_str = c("--CODE_PATH"), 
                          type = "character", 
                          default = CODE_PATH, 
                          metavar = "character", 
                          help = "Path to R Code for Project"),
    
    optparse::make_option(opt_str = c("--FINAL_OUTPUT_PATH"), 
                          type = "character", 
                          default = FINAL_OUTPUT_PATH, 
                          metavar = "character", 
                          help = "Path to Final Output From Code"),

    optparse::make_option(opt_str = c("--LOG_OUTPUT_PATH"), 
                          type = "character", 
                          default = LOG_OUTPUT_PATH, 
                          metavar = "character", 
                          help = "Path to Log Outputs"),

    optparse::make_option(opt_str = c("--SUMMARY_OUTPUT_PATH"), 
                          type = "character", 
                          default = SUMMARY_OUTPUT_PATH, 
                          metavar = "character", 
                          help = "Path to Save the Summary Report")
    
)

if(DO_LOGGING){ 
    log_info("Setting Script Parameters....")
} else {
    message("Setting Script Parameters....")
}

opt_parser <- optparse::OptionParser(option_list = option_list, add_help_option = TRUE)
# opt_parser

Opt_DataPrep <- optparse::parse_args(opt_parser)
# Opt_DataPrep

#################################################################################################
### IMPORT DATA (LATEST FILE) 

if(DO_LOGGING){ 
    log_info("Finding Most Recent File....")
} else {
    message("Finding Most Recent File.....")
}

all_files <- file.info(list.files(Opt_DataPrep$RAW_DATA_PATH, full.names = TRUE))

filename <- rownames(all_files)[which.max(all_files$mtime)]

if(DO_LOGGING){ 
    log_info("Import Data....")
} else {
    message("Import Data....")
}

lyft_dat <- fread(filename)
# lyft_dat

#################################################################################################
### MODIFY DATA

if(DO_LOGGING){ 
    log_info("Modify Data....")
} else {
    message("Import Data....")
}

#colnames(lyft_dat) <- sub(" ", "_", tolower(colnames(lyft_dat)))
# colnames(lyft_dat)
colnames(lyft_dat) <- tolower(gsub("[[:space:]]", "_", colnames(lyft_dat)))
# colnames(lyft_dat)

lyft_dat[, full_name := paste(get(FIRST_NAME_COLUMN_NAME), get(LAST_NAME_COLUMN_NAME), sep = " ")]
# lyft_dat$full_name[1:30]

if(DO_LOGGING){ 
    log_info("The Dataset Includes { nrow(lyft_dat) } Rows")
    log_info("The Dataset Includes Data on { length(unique(lyft_dat$full_name)) } Patients")
} else {
    message("The Dataset Includes Rows: ", nrow(lyft_dat))
    message("The Dataset Includes Data on Patients: ", length(unique(lyft_dat$full_name)))
}

#################################################################################################
### CREATE DATA DICTIONARY

if(DO_LOGGING){ 
    log_info("Create Data Dictionary....")
} else {
    message("Create Data Dictionary....")
}

data_dict <- lyft_dat[, .N, by = .(passenger_number, full_name)][order(-N)]
# data_dict

#################################################################################################
### CREATE DATA TABLE TO STORE SUMMARY REPORT 

SUMMARY_REPORT_OUTPUT <- list() 

#################################################################################################
### CREATE DIRECTORY TO STORE OUTPUS

dir.create(file.path(Opt_DataPrep$FINAL_OUTPUT_PATH, TIMESTAMP_SAVE), showWarnings = FALSE)

NEW_OUTPUT_PATH <- file.path(Opt_DataPrep$FINAL_OUTPUT_PATH, TIMESTAMP_SAVE)

#################################################################################################
### GET DATA FOR EACH MEMBER AND GENERATE OUTPUT  

if(DO_LOGGING){ 
    log_info("Looping Through Each Patient Name....")
} else {
    message("Looping Through Each Patient Name....")
}

# each_name = 2
for(each_name in 1:nrow(data_dict)){

    SUMMARY_REPORT_OUTPUT_TMP <- data.table()  
  
    if(DO_LOGGING){ 
        log_info("Executing: { each_name }")
        log_info("Executing: {data_dict[each_name, full_name]}")
    } else {
        message("Executing: ", each_name)
        message("Executing: ", data_dict[each_name, full_name])
    }
    
    lyft_dat_sub <- lyft_dat[full_name == data_dict[each_name, full_name], ]
    # lyft_dat_sub
  
    patient_name <- unique(lyft_dat_sub[, full_name])

    dir.create(file.path(NEW_OUTPUT_PATH, patient_name), showWarnings = FALSE)
    #dir.create(file.path(MAIN_PATH, "final_output", patient_name), showWarnings = FALSE)
    FINAL_OUTPUT_PATH_PATIENT <- file.path(NEW_OUTPUT_PATH, patient_name)

    SUMMARY_REPORT_OUTPUT_TMP[, file_name := basename(filename)]
    
    SUMMARY_REPORT_OUTPUT_TMP[, passenger_id := data_dict[each_name, get(PASSENGER_NUMBER_COLUMN)]]
    SUMMARY_REPORT_OUTPUT_TMP[, patient_name := patient_name]
    SUMMARY_REPORT_OUTPUT_TMP[, number_of_trips := nrow(lyft_dat_sub)]

    SUMMARY_REPORT_OUTPUT_TMP[, num_unique_days := uniqueN(weekdays(mdy(lyft_dat_sub[, get(PICKUP_DATE_COLUMN)])))]
    SUMMARY_REPORT_OUTPUT_TMP[, num_unique_cities := uniqueN(lyft_dat_sub[, get(PICKUP_CITY_COLUMN)])]
    SUMMARY_REPORT_OUTPUT_TMP[, min_transaction_amount := min(lyft_dat_sub[, get(TRANSACTION_COLUMN)])]
    SUMMARY_REPORT_OUTPUT_TMP[, max_transaction_amount := max(lyft_dat_sub[, get(TRANSACTION_COLUMN)])]
    SUMMARY_REPORT_OUTPUT_TMP[, mean_transaction_amount := mean(lyft_dat_sub[, get(TRANSACTION_COLUMN)], na.rm = TRUE)]
    SUMMARY_REPORT_OUTPUT_TMP[, mean_distance := mean(lyft_dat_sub[, get(DISTANCE_COLUMN)], na.rm = TRUE)]
    SUMMARY_REPORT_OUTPUT_TMP[, mean_duration := mean(lyft_dat_sub[, get(DURATION_COLUMN)], na.rm = TRUE)]
    SUMMARY_REPORT_OUTPUT_TMP[, num_unique_requesters := uniqueN(lyft_dat_sub[, get(REQUESTER_COLUMN)])]

    # output <- lyft_dat_sub[, .(`pickup_date (local)`, `pickup_time (local)`, first_name, 
    #                  last_name, request_address, destination_address, `distance_(miles)`, 
    #                  `duration_(minutes)`, transaction_amount, passenger_number, 
    #                  passenger_number, requester_name)]

    if(DO_LOGGING){ 
        log_info("Executing: Collect Data for Output")
    } else {
        message("Executing: Collect Data for Output")
    }
    
    output <- lyft_dat_sub[, mget(SAVE_THESE_COLUMN_NAMES)]
    # dim(output)             

    if(DO_LOGGING){ 
        log_info("Add Row With Sum Target Variable")
    } else {
        message("Add Row With Sum Target Variable")
    }
    
    sum_target_value = sum(output[, get(TARGET_VARIABLE)])
    
    SUMMARY_REPORT_OUTPUT_TMP[, total_cost := sum_target_value]
    
    output_tmp <- setDT(janitor::adorn_totals(output))
    # output_tmp
    
    #output_tmp[nrow(output_tmp)]
    output_tmp[nrow(output_tmp), setdiff(colnames(output_tmp)[sapply(output_tmp, is.numeric)], 
                                         c(TARGET_VARIABLE, DISTANCE_COLUMN, DURATION_COLUMN)) := ""]
    # output_tmp
    
    output_tmp[, transaction_amount := paste("$", transaction_amount, sep="")]
    
    output_final <- copy(output_tmp)
    
    if(DO_LOGGING){ 
        log_info("The Filtered Dataset Has { nrow(output) } rows")
    } else {
        message("The Filtered Dataset Has rows: ", nrow(output))
    }
    
    if(each_name == 1 && any(SAVE_METHOD %in% "kable")){ 
        if(is.null(webshot:::find_phantom())){
            webshot::install_phantomjs(force = TRUE)
        }
    }
    
    # EACH_SAVE_METHOD = "xls"
    for(EACH_SAVE_METHOD in SAVE_METHOD){

        if(DO_LOGGING){ 
            log_info("Saving Method: { EACH_SAVE_METHOD }")
        } else {
            message("Saving Method: ", EACH_SAVE_METHOD)
        }
      
        tryCatch({
    
            if(EACH_SAVE_METHOD == "kable"){ 
            
                output_df <- as.data.frame(output_final)
                output_html <- kable(output_df, "html")
                output_html_final <- kable_styling(output_html)
                
                save_kable(output_html_final, file.path(NEW_OUTPUT_PATH,
                                                        patient_name,
                                       paste(gsub("[[:space:]]", "_", patient_name), 
                                             TIMESTAMP_SAVE, 
                                             "Final_Output.pdf", sep = "_")))
                
            } else if(EACH_SAVE_METHOD == "pdf"){ 
              
                output_df <- as.data.frame(output_final)
    
                pdf(file.path(NEW_OUTPUT_PATH,
                              patient_name,
                                       paste(gsub("[[:space:]]", "_", patient_name), 
                                             TIMESTAMP_SAVE, 
                                             "Final_Output.pdf", sep = "_")),
                       width = 25, height = 25)
                
                grid.table(output_df)
                
                dev.off()
                
            } else if(EACH_SAVE_METHOD == "excel"){ 
              
                output_dt <- setDT(output_final)
                
                fwrite(output_dt, file.path(NEW_OUTPUT_PATH,
                                            patient_name,
                                       paste(gsub("[[:space:]]", "_", patient_name),
                                             TIMESTAMP_SAVE, 
                                             "Final_Output.xls", sep = "_")))
                
            } else if(EACH_SAVE_METHOD == "csv"){ 
              
                output_dt <- setDT(output_final)
                
                fwrite(output_dt, file.path(NEW_OUTPUT_PATH,
                                            patient_name,
                                       paste(gsub("[[:space:]]", "_", patient_name), 
                                             TIMESTAMP_SAVE, 
                                             "Final_Output.csv", sep = "_")))
                
            }
          
            SUMMARY_REPORT_OUTPUT_TMP[, output_created := "YES"]
            SUMMARY_REPORT_OUTPUT_TMP[, save_method := EACH_SAVE_METHOD]
    
        
        }, error = function(cond) {
    
              SUMMARY_REPORT_OUTPUT_TMP[, output_created := "NO"]
              SUMMARY_REPORT_OUTPUT_TMP[, save_method := EACH_SAVE_METHOD]

              message("eror", cond, sep = "  -  ")
        
        })
    }
  
    if(DO_LOGGING){ 
        log_info("Saving Summary Report Results for { patient_name }")
    } else {
        message("Saving Summary Report Results for: ", patient_name)
    }

    SUMMARY_REPORT_OUTPUT[[each_name]] <- SUMMARY_REPORT_OUTPUT_TMP
    
    Sys.sleep(1)
    
    if(DO_LOGGING){ 
        log_info(" ")
    } else {
        message(" ")
    }

}



if(DO_LOGGING){ 
    log_info("Loop Executed Succesfully")
} else {
    message("Loop Executed Succesfully")
}


SUMMARY_REPORT <- rbindlist(SUMMARY_REPORT_OUTPUT)

SUMMARY_REPORT[, save_method := NULL]

# SUMMARY_REPORT

#################################################################################################
### SAVE SUMMARY RESULTS

if(DO_LOGGING){ 
    log_info("Saving Summary Results")
} else {
    message("Saving Summary Results")
}

if(GENERATE_SUMMARY_OUTPUT){
  
    fwrite(SUMMARY_REPORT, 
           file.path(Opt_DataPrep$SUMMARY_OUTPUT_PATH, 
                     paste("Summary_Report_", TIMESTAMP_SAVE, ".csv", sep = "")))  

}

#################################################################################################
#################################################################################################

end_time <- Sys.time()

time_elapsed <- end_time - start_time

if(DO_LOGGING){ 
    log_appender()
    log_info("Execution Time: { time_elapsed }")
} else {
    message("Script Completed")
    message("Execution Time: ", time_elapsed)
}
 
#################################################################################################
#################################################################################################


