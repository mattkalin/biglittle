# Kalin's Big/Little Pairing Algorithm

Kalin's Big/Little algorithm was made by Matt Kalin in 2020 and was used to make the pairings for Sigma Alpha Mu's Spring 2021 big/little at the University of Maryland.

### Use
First, make sure you have RStudio installed: https://www.rstudio.com/products/rstudio/download/
Then, do the following:
* Save the preferences excel sheet on your computer and update the folder.path and file.name variables in the R file accordingly (line 582)
* Make sure the preferences excel sheet is formatted like the Example data.xlsx file in this repo
  * One sheet titled "Brother" and the other sheet titled "Pledge"
  * Each name must be unique and spelled the same way throughout the sheet
    * In the event of duplicate names, consider using their last name initial (ex: Josh A and Josh M)
  * Headings for the Brother sheet must be: Brother, Pledge 1, Pledge 2, etc
    * Pledge 1 represent's the brother's first choice for a little, Pledge 2 represents his second choice, etc
    * If each brother is given 3 choices (as in the Example data), they will be listed in columns Pledge 1, Pledge 2, and Pledge 3
    * Similarly, headings for the pledge sheet must be: Pledge, Brother 1, Brother 2, etc
* Set who is eligible to take twins, if anyone
  * When a brother gets two littles, that is called getting "twins"
  * This is controlled by the twins.elig variable (line 596)
    * If no one is going to get twins (this should happen when the number of bigs and littles are identical), set twins.elig = c()
    * If every big is eligible to get twins, set twins.elig = "all"
    * If one big is eligible to get twins, specify like twins.elig = "Howard"
    * If multiple specific bigs are eligible to get twins, specify like twins.elig = c("Howard", "George")
* Once you have done all of this, run the algorithm by highlighting from the beginning of the R file through the comment "where the magic happens" (line 698) and clicking run (or command+Enter)



Contact: kalinmatt4@gmail.com


Watch the videos:
1. Introduction
