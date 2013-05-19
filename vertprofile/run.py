#!/usr/bin/env python

import subprocess
import os
import multiprocessing
import socket



cwd = os.getcwd()
def nclprocess(idx,sEdison):
	print "processing model number {}".format(idx)
	with open(os.devnull,'w') as devnull:
		sCmd = "ncl bExecuteOnEdison={} iFileIndex={} {}/vertprof.ncl".format(sEdison,idx,cwd)
		#pid = subprocess.Popen(sCmd,shell=True,stdout=devnull)
		pid = subprocess.Popen(sCmd,shell=True)
		pid.wait()

def run(nModels,sEdison):
	pool = multiprocessing.Pool()
	for i in range(nModels):
		pool.apply_async(nclprocess, args=(i,sEdison))
	pool.close()
	pool.join()

if __name__ == '__main__':
	nModels = input("number of models? ")

	sEdison = 'True' if socket.gethostname().find('edison') != -1 else 'False'
	print "sEdison = {}".format(sEdison)
	run(nModels, sEdison)
	
