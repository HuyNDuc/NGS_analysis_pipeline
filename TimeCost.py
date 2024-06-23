#Created by Huy Duc
# University of Colorado, Anschutz Medical Campus
# FUnctional Genomics Facility

import time as time
import datetime

#function to calculate time
def timeCost(startTime):
    timeCost = time.time() - startTime
    rigtnow = datetime.datetime.now()
    if timeCost < 60:
        timeCost = str("%.2f" % timeCost)
        statement = 'It took '+timeCost+' seconds to finish at %s:%s of %s/%s/%s !' % (rigtnow.hour, rigtnow.minute, rigtnow.day, rigtnow.month, rigtnow.year)
        return statement
    elif timeCost >= 60 and timeCost < 3600:
        minute = timeCost//60
        second = timeCost % 60
        minute = str("%.0f" % minute)
        second = str("%.2f" % second)
        statement = 'It took '+minute + ' minutes, and '+second+' seconds to finish at %s:%s of %s/%s/%s !' % (rigtnow.hour, rigtnow.minute, rigtnow.day, rigtnow.month, rigtnow.year)
        return statement
    elif timeCost >= 3600 and timeCost <86400:
        hour = timeCost//3600
        minute = timeCost % 3600
        minute = minute //60
        second = timeCost % 60
        hour = str("%.0f" % hour)
        minute = str("%.0f" % minute)
        second = str("%.2f" % second)
        statement = 'It took '+hour+' hours, '+ minute + ' minutes and '+second+' seconds to finish at %s:%s:%s of %s/%s/%s !' % (rigtnow.hour, rigtnow.minute, rigtnow.second, rigtnow.day, rigtnow.month, rigtnow.year)
        return statement
    else:
        day = timeCost //86400
        hour = timeCost % 86400 // 3600
        minute = (timeCost % 86400 % 3600)//60
        second = (timeCost % 86400 % 3600 % 60)
        day = str("%.0f" % day)
        hour = str("%.0f" % hour)
        minute = str("%.0f" % minute)
        second = str("%.1f" % second)
        statement = 'It took ' +day + ' days, '+ hour + ' hours, ' + minute + ' minutes and, ' + second + ' seconds to finish at %s:%s of %s/%s/%s !' % (rigtnow.hour, rigtnow.minute, rigtnow.day, rigtnow.month, rigtnow.year)
        return statement


def timestart():
    '''Function to display when the main function starts'''
    rightnow = datetime.datetime.now()
    statement = 'Function started at: %s:%s of %s/%s/%s!' % (rightnow.hour,rightnow.minute,rightnow.day,rightnow.month,rightnow.year)
    return statement