import os
import time
import copy

##########################################################################
class Parallel(object):
    """A class for robustly handling parallel calls to standalone scripts.

    Parameters
    ----------
    max_cpu : int
        Maximum number of parallel threads.
    sleep_interval : int
        Number of seconds to wait between to polls for available threads.
    """

    #------------------------------------------------------------------------#
    def __init__(self, max_cpu, sleep_interval):
        """Constructor of the Parallel class.
        """

        self.max_cpu = max_cpu
        self.sleep_interval = sleep_interval
        self.cmds = []
        self.running = []
    #------------------------------------------------------------------------#

    #------------------------------------------------------------------------#
    def init_cmd_pool(self, cmds):
        """Initializes the command pool.

        Parameters
        ----------
        cmds : List[str]
            List of commands to be issued in parallel.
        
        Returns
        -------
        None
        """

        self.cmds = copy.copy(cmds)
    #------------------------------------------------------------------------#

    #------------------------------------------------------------------------#
    def spawn_ret_pid(self, cmd):
        """Calls the command and returns the process's pid.
        
        Parameters
        ----------
        cmd : str
            The command to run.

        Returns
        -------
        int
            pid of the worker process

        Notes
        -----
        Taken from stackexchange (https://stackoverflow.com/questions/20218570/
        how-to-determine-pid-of-process-started-via-os-system)
        """
    
        print cmd
    
        cmd += " & echo $! > pid"
    
        os.system(cmd)
        pid_file = open("pid", "r")
        pid = int(pid_file.read())
        pid_file.close()
    
        print "pid ", pid
        os.system("rm pid")
    
        return pid
    #------------------------------------------------------------------------#

    #------------------------------------------------------------------------#
    def poll_pids(self, clean = True):
        """Polls the list of running pids and returns the number of them
        that are not alive.

        Parameters
        ----------
        clean : Optional[bool], deafult = True
            Specifies if the dead pids should be popped.

        Returns
        -------
        int
            Number of dead pids.
        """
        
        #........................................................................#
        def check_pid_alive(pid):
            """Checks if the pid is alive.

            Parameters
            ----------
            pid : int
               pid to check.

            Returns
            -------
            bool
                True if alive else False

            Notes
            -----
            Taken from stackexchange (https://stackoverflow.com/questions/568271/
            how-to-check-if-there-exists-a-process-with-a-given-pid-in-python)
            """

            try:
                os.kill(pid, 0)
            except OSError:
                return False
            else:
                return True
        #........................................................................#

        removal_list = []
        for pid in self.running:
            if not check_pid_alive(pid):
                removal_list.append(pid)
    
        if clean:
            for pid in removal_list:
                self.running.remove(pid)

        return len(removal_list)            
    #------------------------------------------------------------------------#
    
    #------------------------------------------------------------------------#
    def run(self):
        """Runs the initialized pool.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        
        i = 0
        done = False
        while i < len(self.cmds) or done:
            if done:
                while self.running:
                    self.poll_pids()
                break
            while len(self.running) < self.max_cpu:
                self.running.append(self.spawn_ret_pid(self.cmds[i]))
                i += 1
                if i == len(self.cmds):
                    done = True
                    break
            if done:
                continue
            freed = 0
            while freed == 0:
                freed = self.poll_pids()
                time.sleep(self.sleep_interval)
    #------------------------------------------------------------------------#
##########################################################################
