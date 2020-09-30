#progress bar class, has numpy, time and datetime dependencies
import numpy as np
import time
import datetime

class ProgressBar(object):
    """
    Call in a loop to create terminal progress bar
    """
    
    def __init__(self, total, decimals = 1, length = 50, fill = '=', tip = '>', end_string = '\r', 
                 comp_line = '', global_time_measure = False, iteration_time_measure = False, flush_value = False):
        """
        Initializes the progress bar
        
        Inputs:
        total                  - Required  : total iterations (Int)
        decimals               - Optional  : positive number of decimals in percent complete (Int)
        length                 - Optional  : character length of bar (Int)
        fill                   - Optional  : bar fill character (Str)
        tip                    - Optional  : bar tip character (Str)
        end_string             - Optional  : end character (e.g. "\r", "\r\n") (Str)
        comp_line              - Optional  : line to print after completion (Str)
        global_time_measure    - Optional  : determines if the total runtime is measured (Bool)
        iteration_time_measure - Optional  : determines if the runtime of a single loop is measured (Bool)
        flush_value            - Optional  : flush argument inside each print command (Bool)
        """
        
        self.total = total
        
        self.decimals = decimals
        self.length = length
        self.fill = fill
        self.tip = tip
        self.end_string = end_string 
        self.comp_line = comp_line
        self.global_time_measure = global_time_measure
        self.iteration_time_measure = iteration_time_measure
        self.flush_value = flush_value
            
        if self.iteration_time_measure == True:
            self.iteration_time_list = []
            self.start_iteration_time = time.time()
        
        if self.global_time_measure == True:
            self.start_global_time = time.time()
                
    def start_iteration_measure(self):
        """
        Starts counting iteration time
        """
        if self.iteration_time_measure != True:
            raise ValueError('The iteration time measure has to be initialized!')
        self.start_iteration_time = time.time()
                                                                         
    def show(self, iteration, prefix='', suffix=''):
        """
        Prints the progress bar
        
        Inputs:
        iteration   - Required  : current iteration (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)         
        """
        
        self.prefix = prefix
        self.suffix = suffix
        
        percent = ("{0:." + str(self.decimals) + "f}").format(100 * (iteration / float(self.total)))
        filled_length = int(self.length * iteration // self.total)
        
        if iteration < self.total:
            if self.iteration_time_measure == True:
                iteration_time = time.time() - self.start_iteration_time
                self.iteration_time_list.append(iteration_time)
                estimate_time = np.round(np.mean(self.iteration_time_list)*(self.total-iteration), decimals=0)
                self.suffix = '- {} left'.format(datetime.timedelta(seconds=estimate_time)) + suffix 
                
            bar = self.fill * (filled_length - 1) + self.tip + '.' * (self.length - filled_length)
            print(f'\r{self.prefix} [{bar}] {percent}% {self.suffix}', end = self.end_string, flush = self.flush_value)
            
        # Print New Line on Complete and full bar without tip
        if iteration == self.total:
            if self.global_time_measure == True:
                seconds_global_time = np.round(time.time() - self.start_global_time, decimals=0)
                self.global_time = datetime.timedelta(seconds=seconds_global_time)
                self.suffix = '- {} elapsed'.format(self.global_time) + suffix
                
            bar = self.fill * filled_length + '.' * (self.length - filled_length)
            print(f'\r{self.prefix} [{bar}] {percent}% {self.suffix}', end = '', flush = self.flush_value)
            print(self.comp_line)