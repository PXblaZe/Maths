from curses import *
import time

def test(stdscr):
    stdscr.addstr(0, 0, 'Hello World')
    stdscr.refresh()
    time.sleep(3)

wrapper(test)