require(data.table)

source("bin/All_Classes.R")
source("bin/obsolete.R")
source("bin/ora.R")
source("bin/taxonomy.R")

# Upload dataset ----

convertColumns <- c("UNIPROTID", "ENTREZID", "SYMBOL", "GENEID", "PROTEINID")
assign("convertColumns", convertColumns, envir = .GlobalEnv)


# User timeout ----

timeoutMinutes <- 20

inactivity <- sprintf("function idleTimer() {
var t = setTimeout(logout, %s);
window.onmousemove = resetTimer; // catches mouse movements
window.onmousedown = resetTimer; // catches mouse movements
window.onclick = resetTimer;     // catches mouse clicks
window.onscroll = resetTimer;    // catches scrolling
window.onkeypress = resetTimer;  //catches keyboard actions

function logout() {
Shiny.setInputValue('timeOut', '%s minutes of')
}

function resetTimer() {
clearTimeout(t);
t = setTimeout(logout, %s);  // time is in milliseconds (1000 is 1 second)
}
}
idleTimer();", timeoutMinutes*6e4, timeoutMinutes, timeoutMinutes*6e4)
