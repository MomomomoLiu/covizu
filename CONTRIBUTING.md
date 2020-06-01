#The following Environment Variables need to be defined for GISAID downloads

gisaid_pw_variable=
gisaid_u_variable=

#Downloads can be automated through crontab on linux:
#Bash environment needs to be initialized in this process

0 0 * * * nohup python /home/covizu/SeliniumAutobot.py >> /home/covizu/Autobot.log 2>&1
