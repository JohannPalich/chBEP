#!/usr/bin/env python
# -*- coding: utf-8 -*-

#usual psychopy import
from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import visual, core, data, event, logging, gui, monitors, misc, filters
from psychopy.constants import *  # things like STARTED, FINISHED
from psychopy.tools.monitorunittools import deg2pix, pix2deg
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
from random import choice, uniform, sample
#specific imports for this script
import codecs, math
from math import cos, sin
import csv
import copy
#import analysis

from heapq import nsmallest
from operator import itemgetter

# Store info about the experiment session
expName = u'cb_v2'  # from the Builder filename that created this script
expInfo = {u'1. Ваше имя (или ник)': u'', u'2. Возраст': u'', u'3. Пол (m/f)': u'', u'4. Доминирующий глаз (l/r)':u'',u'5. Метод трекинга (l/r/b)':u'', 'Tracking Enabled':'0'}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
tracking_mode = expInfo[u'5. Метод трекинга (l/r/b)']
tracking_en = int(expInfo['Tracking Enabled'])


if dlg.OK == False or tracking_mode not in ['l','r','b']: core.quit()    # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName


def grabScreenshot(file, win):
    win.getMovieFrame()
    filename = str(file) + '.png'
    win.saveMovieFrames(filename)
    win.movieFrames = []

# ---------------------------------------------
#---- connect to iView
# ---------------------------------------------

if tracking_en:
    from iViewXAPI import *
    res = iViewXAPI.iV_SetLogger(c_int(1), c_char_p("iViewXSDK_Python_GazeContingent_Demo.txt"))
    res = iViewXAPI.  iV_Connect(c_char_p('127.0.0.1'), c_int(4444), c_char_p('127.0.0.1'), c_int(5555))

    res = iViewXAPI.iV_GetSystemInfo(byref(systemData))
    print "iV_GetSystemInfo: " + str(res)
    print "Samplerate: " + str(systemData.samplerate)
    print "iViewX Verion: " + str(systemData.iV_MajorVersion) + "." + str(systemData.iV_MinorVersion) + "." + str(systemData.iV_Buildnumber)
    print "iViewX API Verion: " + str(systemData.API_MajorVersion) + "." + str(systemData.API_MinorVersion) + "." + str(systemData.API_Buildnumber)


# ---------------------------------------------
#---- configure and start calibration
# ---------------------------------------------

def calibrate_and_validate():
    displayDevice = 0
    dev_x = 5
    dev_y = 5
    while dev_x > 0.5 or dev_y > 0.5:
        calibrationData = CCalibration(5, 1, displayDevice, 0, 1, 20, 239, 1, 10, b"")

        res = iViewXAPI.iV_SetupCalibration(byref(calibrationData))
        print "iV_SetupCalibration " + str(res)
        res = iViewXAPI.iV_Calibrate()
        print "iV_Calibrate " + str(res)

        res = iViewXAPI.iV_Validate()
        print "iV_Validate " + str(res)
        res = iViewXAPI.iV_GetAccuracy(byref(accuracyData), 0)
        print "iV_GetAccuracy " + str(res)
        print "deviationXLeft " + str(accuracyData.deviationLX) + " deviationYLeft " + str(accuracyData.deviationLY)
        print "deviationXRight " + str(accuracyData.deviationRX) + " deviationYRight " + str(accuracyData.deviationRY)
        dev_x = float(accuracyData.deviationLX) if tracking_mode=='l' or tracking_mode=='b' else float(accuracyData.deviationRX)
        dev_y = float(accuracyData.deviationLY) if tracking_mode=='l' or tracking_mode=='b' else float(accuracyData.deviationRY)
    
if tracking_en:
    calibrate_and_validate()

globalClock = core.Clock()

# Setup files for saving
if not os.path.isdir('data'):
    os.makedirs('data')  # if this fails (e.g. permissions) we will get error

filename = 'data' + os.path.sep + '%s' % (expInfo['date'])
logFile = logging.LogFile(filename + '.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
                                 extraInfo=expInfo, runtimeInfo=None,
                                 originPath=None,
                                 savePickle=True, saveWideText=True,
                                 dataFileName=filename)

# Setup the Window
monitorName = 'LG Flatron L1718S'
monitorName = 'SyncMaster 797MB'
monitorName = 'ASUS ML 239'
#monitorName = 'HP Compaq LA2205wg'
monitorName = 'Dell Vostro 5470'
monitorName = 'testMonitor'
frameRate = 60
mon = monitors.Monitor(monitorName)

win = visual.Window(size=[1024, 768], fullscr=True, screen=0, allowGUI=True, allowStencil=True,
                    monitor=mon, color=[1, 1, 1], colorSpace=u'rgb', units=u'pix')

screen_dim_deg = [pix2deg(win.size[0], mon), pix2deg(win.size[1], mon)]
#print 'Refresh rate %.3f' % win.getActualFrameRate()

f = codecs.open("instr.txt", "r", encoding='utf-8')

try:
    f = codecs.open("instr.txt", "r", encoding='utf-8')
    try:
        instrContent = [f.read()]
    finally:
        f.close()
except IOError:
    pass

myMouse = event.Mouse(win=win)
instr_text = visual.TextStim(win=win, name='instr_text',
                             text=instrContent[0], font=u'Arial', color='black',
                             units=u'norm', pos=[0, 0], height=0.05, wrapWidth=1.8)

instr_text.draw()
win.flip()
myMouse.clickReset()
buttons = [0]

while not buttons[0]:
    buttons = myMouse.getPressed()
    core.wait(0.05)


fp=visual.Polygon(win=win, name='polygon',
        edges = 96, size=[1,1],
        ori=45, pos=[0,0],
        lineWidth=0, lineColor='red', lineColorSpace=u'rgb',
        fillColor='red', fillColorSpace=u'rgb',
        opacity=1,interpolate=True,units='cm')

progress_bar_frame = visual.Rect(win, width=1, height=0.02, units=u'norm', pos=[0, -0.95], fillColor=None,
                                 lineColor='#aaaaaa', autoLog=False)
progress_bar = visual.Rect(win, width=1, height=0.02, units='norm', pos=[0, -0.95], fillColor='#aaaaaa',
                           lineColor='#aaaaaa', autoLog=False)

matSize=(5,5)
imSize=2.5/1.5 #70
imSize_70 = imSize
imSize_94 = 94*imSize/70 #image size for rotated images is larger
nImages=matSize[0]*matSize[1]

myClock=core.Clock()
probeClock = core.Clock()
print 'nImages %g' %nImages

posList=[]
for x in xrange(matSize[0]):
    for y in xrange(matSize[1]):
        posList.append(((x-matSize[0]/2+1/2)*imSize*1.5, (y-matSize[1]/2+1/2)*imSize*1.5))

imList= []
imList_changed=[]

def gen_stimuli():
    colors_list=['blue','brown']
    rot_list=['0','45']

    shapes_list=range(1,65)

    shuffle(colors_list)
    shuffle(rot_list)

    colors=colors_list*(int(nImages/len(colors_list))+1)
    rots=rot_list*(int(nImages/len(rot_list))+1)

    shuffle(colors)
    shuffle(rots)
    shuffle(shapes_list)
#    print colors
    parMatrix=[]
    for x in xrange(matSize[0]):
        for y in xrange(matSize[1]):
#            print (x, y)
            parMatrix.append([colors.pop(), rots.pop(),shapes_list[x*5+y], posList[x*5+y]])

    change_type=trial['change_type']
    change_pos=trial['change_pos']
    if len(imList)>0:
        for i in xrange(len(parMatrix)):
            imList[i].image = '%s_%s/%03d.png' % tuple(parMatrix[i][0:3])
            imList[i].size = imSize_70 if parMatrix[i][1]=='0' else imSize_94
            imList_changed[i].image = '%s_%s/%03d.png' % tuple(parMatrix[i][0:3])
            imList_changed[i].size = imSize_70 if (parMatrix[i][1]=='0' and (change_pos!=i or change_type!='rot')) or (parMatrix[i][1]!='0' and change_type=='rot' and change_pos==i) else imSize_94
    else:
        for i in xrange(len(parMatrix)):
            imList.append(visual.ImageStim(win, '%s_%s/%03d.png' % tuple(parMatrix[i][0:3]),
                                                 pos=(parMatrix[i][3]),
                                                 units='deg', size=imSize_70 if parMatrix[i][1]=='0' else imSize_94))
            imList_changed.append(visual.ImageStim(win, '%s_%s/%03d.png' % tuple(parMatrix[i][0:3]),
                                                 pos=(parMatrix[i][3]),
                                                 units='deg', size=imSize_70 if (parMatrix[i][1]=='0' and (change_pos!=i or change_type!='rot')) or (parMatrix[i][1]!='0' and change_type=='rot' and change_pos==i) else imSize_94))
        
#    screenshot.draw()
    #win.flip()
    #screenshot = visual.BufferImageStim(win, stim=imList)

#    print parMatrix[change_pos]
    #print imList[change_pos].image
    with open(filename+'_par_trial_'+`trials.thisN`+'.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerows([['color','rot','shape','pos']])
        writer.writerows(parMatrix)

    cp_orig=parMatrix[change_pos]
    trials.addData('target_orig_color',cp_orig[0])
    trials.addData('target_orig_rot',cp_orig[1])
    trials.addData('target_orig_shape',cp_orig[2])
    if change_type=='color':
        if parMatrix[change_pos][0]=='blue':
            parMatrix[change_pos][0]='brown'
        else:
            parMatrix[change_pos][0]='blue'
    elif change_type=='rot':
        if parMatrix[change_pos][1]=='0':
            parMatrix[change_pos][1]='45'
        else:
            parMatrix[change_pos][1]='0'
    elif change_type=='shape':
        parMatrix[change_pos][2]=shapes_list[nImages+1]

    cp_new=parMatrix[change_pos]
    trials.addData('target_orig_color',cp_new[0])
    trials.addData('target_orig_rot',cp_new[1])
    trials.addData('target_orig_shape',cp_new[2])

#    print 'Change Position %i'%change_pos

#    print 'Change Type %s'%change_type
#    print parMatrix[change_pos]
    #print imList[change_pos].image

    imList_changed[change_pos].setImage('%s_%s/%03d.png' % tuple(parMatrix[change_pos][0:3]))
    #print imList[change_pos]._imName
    #screenshot_changed = visual.BufferImageStim(win, stim=imList)

    return (parMatrix)



probe = visual.Circle(win, 0.2, 96, units='deg',  fillColor='#a8a8a8', lineWidth=0)

    
maskTime=5
showTime=15
fpTime=1 #in second
change_types=['color','rot','shape']
trainingList = data.createFactorialTrialList(factors={'change_type': change_types,'change_pos':range(nImages)})
shuffle(trainingList)

fullPosList=range(nImages)
shuffle(fullPosList)

mainTrials=[]
ecc_list=[x - .5 for x in range(1, 8)]
angle_list=range(0,360,30)


#for i in change_types:
#    print i
#    print sample(range(nImages),5)
factorsList = data.createFactorialTrialList(factors={'change_type': change_types, 'trial_probe_ecc':ecc_list, 'trial_probe_angle':angle_list})
    #print 'factorsList: %g'%len(factorsList)
    #mainTrials=mainTrials+factorsList

shuffle(mainTrials)

trials = data.TrialHandler(nReps=1.0, method='random', extraInfo=expInfo, trialList=factorsList)
thisExp.addLoop(trials)

print 'Total trials: %g'%trials.nTotal


# aperture=visual.Rect(win=win,width=1, height=1, fillColor='red', units='norm')
# gaussTexture = filters.makeMask(196, shape='gauss')

# invGaussTexture = -gaussTexture #make opaq->transp and vice versa

def replicate(a, xy, se, n):
    rptIdx = np.ones(a.shape[0 if xy == 'X' else 1], dtype=int)
    rptIdx[0 if se == 'start' else -1] = n + 1
    return np.repeat(a, rptIdx, axis=0 if xy == 'X' else 1)

#print np.shape(invGaussTexture)

# invGaussTexture = replicate(invGaussTexture,'X','start',30)
# invGaussTexture = replicate(invGaussTexture,'X','end',30)
# invGaussTexture = replicate(invGaussTexture,'Y','start',30)
# invGaussTexture = replicate(invGaussTexture,'Y','end',30)

#print np.shape(invGaussTexture)

helper = visual.TextStim(win=win, text='', units=u'norm', pos=[0, 0.9], height=0.05, autoLog=False, color='black')
helper2 = visual.TextStim(win=win, text=u'Если вы заметили меняющуюся картинку, и она сейчас показана, \nщелкните на нее правой кнопкой.', units=u'norm', pos=[0, -0.9], height=0.04,alignHoriz='center', autoLog=False, color=[0,0,0])


origMaskSize=256*4
#gaussMask = visual.GratingStim(win,mask=invGaussTexture,tex=None, contrast=1,  size=origMaskSize, units='pix')


#maskWidth=gaussMask.size[0]
mousePosText=visual.TextStim(win, '',pos=[0.9,0.9],height=0.05, units='norm',color='black')

if tracking_en:
    iViewXAPI.iV_StopRecording()  
    res = iViewXAPI.iV_ClearRecordingBuffer()
    print 'Clear buffer %s' % res
#iViewXAPI.iV_SetEventDetectionParameter(200, 30)
prev_calib_time = globalClock.getTime()

for trial in trials:
    if (globalClock.getTime()-prev_calib_time)>(60*15):
        instr_text.text=u'Вы можете передохнуть. Нажмите на клавишу "S" для продолжения.'
        instr_text.draw()
        win.flip() 
        myMouse.clickReset()
        buttons = [0]

#       while not buttons[0]:
#           buttons = myMouse.getPressed()
#           core.wait(0.05)

        while not len(event.getKeys(['s'])):
            core.wait(0.05)  

        win.winHandle.minimize() # minimise the PsychoPy window
        win.fullscr = False # disable fullscreen
        win.flip() # redraw the (minimised) window

        calibrate_and_validate()

        win.winHandle.maximize()
        win.fullscr = True 
        win.winHandle.activate()
        win.flip()
        prev_calib_time=globalClock.getTime()
        
    mousePosText.autoDraw=False
    probes_loop = data.TrialHandler(nReps=999, method='sequential', seed=None, trialList=None, extraInfo=expInfo, name='probes_loop')
    
    thisExp.addLoop(probes_loop)
    frameN = totalFrameN = 0
    while True:
        
        change_pos = choice(range(nImages))
        target_pos = posList[change_pos]
        
        angle_deg = trial['trial_probe_angle']
        angle = np.pi*angle_deg/180
        eccentr = trial['trial_probe_ecc']
        max_probe_y = target_pos[1]+imSize+eccentr*sin(angle) #we add imSize as it is the minimal distance for gaze position to be "on target"
        min_probe_y = target_pos[1]-imSize+eccentr*sin(angle)
        if abs(max_probe_y)<screen_dim_deg[1]/2 and abs(min_probe_y)<screen_dim_deg[1]/2:
            break       
            

    trials.addData('change_pos',change_pos)
    trial['change_pos']=change_pos
    
    parMatrix = gen_stimuli()
    
    
    #print posList[change_pos]#[0]
    #print imList[change_pos].pos
    
    fp.draw()
    progress_bar.setWidth(trials.thisN / trials.nTotal)
    progress_bar_frame.draw()
    progress_bar.draw()
    #screenshot.draw()
#    for img in imList:
#        img.draw()
    win.flip()

    core.wait(fpTime)
    continueRoutine=True
    mouse_positions=[]
    engaged=False
    engaged = 0
    prev_probe_frame = 0
    next_probe = 0
    #sampleData = []
    while next_probe<2 or next_probe>14:
        next_probe = np.random.gamma(8, 0.8, 1)[0]
    
    nPasses=0
    helper.setText(u'Какая картинка меняется?')
    helper.setAutoDraw(True)
    myClock.reset()
    probeClock.reset()
    prev_time=myClock.getTime()
    if tracking_en:
        iViewXAPI.iV_StopRecording()
        iViewXAPI.iV_ClearRecordingBuffer()
    buttons=myMouse.getPressed()
    buttons=[]
    print('Here1')
    if tracking_en:
        iViewXAPI.iV_StartRecording()
    eye_data=[]
    eye_positions=[]
    last_known_eye_pos=[0,0]
    lastEyeTimeStamp = -1
    if tracking_en:
        res = iViewXAPI.iV_GetSample(byref(sampleData))
        eye_start_time=sampleData.timestamp
    first_probe_on_target = 1
    pursuitStartFrame = 0
    myMouse.clickReset() 
    current=-1
    on_target = current_item = item_change_frameN = item_change_totalFrameN = 0 
    
    while continueRoutine:
        mx, my = myMouse.getPos()
        mx=pix2deg(mx,mon)
        my=pix2deg(my,mon)
        mouse_positions.append([mx, my])
        if tracking_en:
            res = iViewXAPI.iV_GetSample(byref(sampleData))
        #res1 = iViewXAPI.iV_GetEvent(byref(eventData))
            if res == 1:
                ly = sampleData.leftEye
                ry = sampleData.rightEye
                lastEyeTimeStamp = sampleData.timestamp
                ly.gazeX = pix2deg(sampleData.leftEye.gazeX - win.size[0]/2, mon)
                ly.gazeY = pix2deg(-1 * (sampleData.leftEye.gazeY - win.size[1]/2), mon)
                ry.gazeX = pix2deg(sampleData.rightEye.gazeX - win.size[0]/2, mon)
                ry.gazeY = pix2deg(-1 * (sampleData.rightEye.gazeY - win.size[1]/2), mon)
                eye_data.append([sampleData.timestamp,eye_start_time, sampleData.timestamp-eye_start_time, 
                round(ly.gazeX, 4), round(ly.gazeY,4), ly.eyePositionZ, ly.diam, round(ry.gazeX,4), round(ry.gazeY,4), 
                ry.eyePositionZ, ry.diam]) 
                eye_pos = [ly.gazeX, ly.gazeY] if ly.diam > 2 and ly.diam < 7 else [ry.gazeX, ry.gazeY] if ry.diam > 2 and ry.diam < 7 else []
                if len(eye_pos)>1:
                    last_known_eye_pos=eye_pos
            else:
                eye_pos=[]
        else:
            eye_pos=[mx, my]
        eye_positions.append(eye_pos)   
        totalFrameN+=1
        
        if len(eye_positions)>40:
            prev_item = current_item
            min_dist = np.linalg.norm(np.array(last_known_eye_pos) - np.array(imList[current_item].pos))
                
            for i in range(nImages):
                dist = np.linalg.norm(np.array(last_known_eye_pos) - np.array(imList[i].pos))
                if dist < min_dist:
                    current_item = i
                    min_dist = dist
            
            if min_dist>imSize: # eye outside the stimuli field
                current_item = -1
            
            if current_item!=prev_item:
                item_change_frameN = frameN
                item_change_totalFrameN = totalFrameN
                pursuitStartFrame = frameN
                probe.autoDraw = False # probe disappears whenever we change the viewed item
            
            on_target = 1 if current_item == change_pos else 0
            
        min_n = 5 # minimal number of frames for the image considered to be seen
        
        # when pursuitStartFrame falls within the last min_n frames of either phase or during mask period, it is not considered as "seeing the stimuli"
        # thus, we need to update it if observer continues to look at the target
        
        if pursuitStartFrame!=0 and (pursuitStartFrame>(showTime-min_n) and pursuitStartFrame<(showTime+maskTime+min_n) and frameN>(showTime+maskTime+min_n)) or (pursuitStartFrame>(2*showTime+maskTime-min_n) and frameN>=min_n and frameN<showTime):
            pursuitStartFrame = frameN
            print('MyLittleProby')
        print(current)

        if current_item!=-1 and ((pursuitStartFrame>=0 and pursuitStartFrame<(showTime-min_n) and frameN>(showTime+maskTime+min_n)) or #eye comes to target during the 1st phase
            (frameN>=min_n and frameN<showTime and pursuitStartFrame>(showTime+maskTime+min_n) and pursuitStartFrame<(2*showTime+maskTime-min_n))): #eye comes to target during the 2nd phase
            seen_both_versions_for_n_frames = 1
        else:
            seen_both_versions_for_n_frames = 0
        
        if (on_target and first_probe_on_target or totalFrameN>(prev_probe_frame+next_probe*frameRate)) and seen_both_versions_for_n_frames and abs(last_known_eye_pos[0]) < imSize*3 and abs(last_known_eye_pos[1]) < imSize*3:
            thisExp.nextEntry()
            trials.addData('pursuitStart',pursuitStartFrame)
            trials.addData('item_change_frameN',item_change_frameN)
            trials.addData('item_change_totalFrameN',item_change_totalFrameN)
            trials.addData('pursuitEnd',frameN)
            trials.addData('pursuitDuration',totalFrameN-item_change_totalFrameN)
            nPasses+=1
            engaged=0
            #mousePosText.setText(nPasses)
#            print('Yay!')
#            print(last_known_eye_pos)
#            print(len(eye_positions))
            probes_loop.addData('probe_over_target',1)
            probes_loop.addData('probe_start_frameN',totalFrameN)
            probes_loop.addData('probe_delay',next_probe)
            probes_loop.addData('probe_start_time',myClock.getTime())
            probes_loop.addData('probe_start_rel_time',probeClock.getTime())
            probes_loop.addData('eyeStartTimeStamp',lastEyeTimeStamp)
            probes_loop.addData('eyeStartX',last_known_eye_pos[0])
            probes_loop.addData('eyeStartY',last_known_eye_pos[1])
            
            prev_probe_frame = totalFrameN
            if first_probe_on_target:
                angle_deg = trial['trial_probe_angle']
                eccentr=trial['trial_probe_ecc']
                angle=np.pi*angle_deg/180
                probe_xy = [last_known_eye_pos[0]+eccentr*cos(angle), last_known_eye_pos[1]+eccentr*sin(angle)]
            else:
                while True:
                    angle_deg = choice(angle_list)
                    angle=np.pi*angle_deg/180
                    eccentr=choice(ecc_list)
                    probe_xy = [last_known_eye_pos[0]+eccentr*cos(angle), last_known_eye_pos[1]+eccentr*sin(angle)]
                    if abs(probe_xy[0])<screen_dim_deg[0]/2 and abs(probe_xy[1])<screen_dim_deg[1]/2:
                        break

            probe.pos = probe_xy
            probes_loop.addData('probe_x',probe_xy[0])
            probes_loop.addData('probe_y',probe_xy[1])
            probes_loop.addData('probe_angle',angle_deg)
            probes_loop.addData('probe_ecc',eccentr)
            probe.autoDraw = True
            #probe.autoDraw = False
            next_probe = 0
            first_probe_on_target = 0 
            while next_probe<2 or next_probe>14:
                next_probe = np.random.gamma(8, 0.8, 1)[0]
            probeClock.reset()
        
        if current_item!=-1 and totalFrameN>(prev_probe_frame+next_probe*frameRate) and seen_both_versions_for_n_frames  and abs(last_known_eye_pos[0]) < imSize*3 and abs(last_known_eye_pos[1]) < imSize*3:
            thisExp.nextEntry()
            trials.addData('pursuitStart',pursuitStartFrame)
            trials.addData('item_change_frameN',item_change_frameN)
            trials.addData('item_change_totalFrameN',item_change_totalFrameN)
            trials.addData('pursuitEnd',frameN)
            trials.addData('pursuitDuration',totalFrameN-item_change_totalFrameN)
            
            #print(last_known_eye_pos)
            #print(len(eye_positions))
            probes_loop.addData('probe_over_target',0)
            probes_loop.addData('probe_start_frameN',totalFrameN)
            probes_loop.addData('probe_delay',next_probe)
            probes_loop.addData('probe_start_time',myClock.getTime())
            probes_loop.addData('probe_start_rel_time',probeClock.getTime())
            probes_loop.addData('eyeStartTimeStamp',lastEyeTimeStamp)
            probes_loop.addData('eyeStartX',last_known_eye_pos[0])
            probes_loop.addData('eyeStartY',last_known_eye_pos[1])
            prev_probe_frame = totalFrameN
            
            while True:
                angle_deg = choice(angle_list)
                angle=np.pi*angle_deg/180
                eccentr=choice(ecc_list)
                probe_xy = [last_known_eye_pos[0]+eccentr*cos(angle), last_known_eye_pos[1]+eccentr*sin(angle)]
                if abs(probe_xy[0])<screen_dim_deg[0]/2 and abs(probe_xy[1])<screen_dim_deg[1]/2:
                    break
            
            probe.pos = probe_xy
            probes_loop.addData('probe_x',probe_xy[0])
            probes_loop.addData('probe_y',probe_xy[1])
            probes_loop.addData('probe_angle',angle_deg)
            probes_loop.addData('probe_ecc',eccentr)
            
            probe.autoDraw = True
            next_probe = 0
            while next_probe<2 or next_probe>14:
                next_probe = np.random.gamma(8, 0.8, 1)[0]
            probeClock.reset()
        
    
#        mousePosText.setText('%i, %i'%(mx,my))
#        mousePosText.color='black'
        
            
            #mousePosText.autoDraw=False
            
            

        # if frameN % 3==0:
            # print win.size
            # print [mx, my]/win.size*2
            # gaussMask.setPos([mx,my])

        if frameN>=0 and frameN<showTime:
            for im in imList:
                pos=im.pos
                dist=np.linalg.norm((mx, my)- pos)
                #if dist<(gaussMask.size[0]/2-128):
                im.draw()
#               im.draw()
                if frameN==0:
                    time=myClock.getTime()
                #logging.debug('A start. Time %.3f, delta: %.3f' %(time, time-prev_time))
                prev_time=time
        elif frameN==(showTime):
            time=myClock.getTime()
            #logging.debug('A stop. Time %.3f, delta: %.3f' %(time, time-prev_time))
            prev_time=time
        elif frameN>=(showTime+maskTime) and frameN<35:
            if frameN==(showTime+maskTime):
                time=myClock.getTime()
                #logging.debug('A\' start. Time %.3f, delta: %.3f' %(time, time-prev_time))
                prev_time=time
            for im in imList_changed:
                pos=im.pos
                dist=np.linalg.norm((mx, my)- pos)
#                print 'dist: %.0f' % dist
#                print 'mouse: %i, %i' % (mx, my)
#                print pos

                #if dist<(gaussMask.size[0]/2-128):
                im.draw()
        elif frameN==35:
            #screenshot_changed.autoDraw=False
            time=myClock.getTime()
            #logging.debug('A\' stop. Time %.3f, delta: %.3f' %(time, time-prev_time))
            prev_time=time
        elif frameN==39:
            time=myClock.getTime()
            #logging.debug('Time %.3f, delta: %.3f' %(time, time-prev_time))
            prev_time=time
            frameN=-1
        #gaussMask.draw()
        #fp.draw()
        progress_bar_frame.draw()
        progress_bar.draw()
#        mousePosText.draw()
        win.flip()
        frameN+=1

        if len(event.getKeys('space')):
            probes_loop.addData('probe_rt_trialstart',myClock.getTime())
            probes_loop.addData('probe_rt_rel',probeClock.getTime())
            probes_loop.addData('probe_rt_frameN', totalFrameN)
            probes_loop.addData('eyeEndTimeStamp',lastEyeTimeStamp)
            probes_loop.addData('eyeEndX',last_known_eye_pos[0])
            probes_loop.addData('eyeEndY',last_known_eye_pos[1])
            probe.autoDraw = False

 
        buttons, times = myMouse.getPressed(True)
        if buttons[0]:
            mx, my = myMouse.getPos()
            mx=pix2deg(mx, mon)
            my=pix2deg(my, mon)
            for i in range(nImages):
                rt=myClock.getTime()
                if imList[i].contains(mx,my):
                    if tracking_en:
                        iViewXAPI.iV_StopRecording()
                    print(imList[i].pos[0])
                    trials.addData('chosen_x',imList[i].pos[0])
                    trials.addData('chosen_y',imList[i].pos[1])
                    trials.addData('mx',mx)
                    trials.addData('my',my)
                    mouse_positions.append([mx, my])
                    trials.addData('answer_rt',rt)
                    trials.addData('answer_frameN', totalFrameN)

                    trials.addData('chosen_pos',i)
                    trials.addData('chosen_color',parMatrix[i][0])
                    trials.addData('chosen_rot',parMatrix[i][1])
                    trials.addData('chosen_shape',parMatrix[i][2])
                    trials.addData('chosen_eye_TimeStamp',lastEyeTimeStamp)
                    correct = 1 if trials.thisTrial['change_pos']==i else 0
                    trials.addData('chosen_correct',correct)
                    trials.addData('chosen_rt',times[0])
                    trials.addData('nPasses',nPasses)
                    continueRoutine=False
                    probe.autoDraw = False
                    win.flip()
                    
                    outputfile = os.getcwd()+'\\'+filename+'_trial_'+`trials.thisN`+'_eye'
                    if tracking_en:
                        res = iViewXAPI.iV_SaveData(str(outputfile),'','', 1)
                        print 'iV_SaveData ' + str(res)
                    print "data saved to: " + outputfile

                    with open(filename+'_trial_'+`trials.thisN`+'.csv', 'wb') as f:
                        writer = csv.writer(f)
                        writer.writerows(mouse_positions)
                    with open(filename+'_trial_'+`trials.thisN`+'_eye.csv', 'wb') as f:
                        writer = csv.writer(f)
                        writer.writerow(['timestamp','trialstarttime','timeafter','lx','ly','lz','ld','rx','ry','rz','rd', 'mx','my','maskSize'])
                        writer.writerows(eye_data)
                        #writer.writerows(eye_positions)

                    with open(filename+'_trial_'+`trials.thisN`+'.csv', 'wb') as f:
                        writer = csv.writer(f)
                        writer.writerows(mouse_positions)
                    mouse_positions=[]

                    break
            myMouse.clickReset()
            buttons=[]
            while myMouse.getPressed()[0]==1:
                pass

        if len(event.getKeys(['escape'])):
            core.quit()
        elif len(event.getKeys(['s'])):
            grabScreenshot('sshot_%i_%i' % (trials.thisN,len(mouse_positions)), win)
        if len(event.getKeys('escape')):
            core.quit()
    thisExp.nextEntry()


win.close()
del thisExp

#analysis.main(expInfo['date'])
