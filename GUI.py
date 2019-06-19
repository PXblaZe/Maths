from curses import *
from curses import textpad
def chrxp(stdscr, ):
    curs_set(0)
    init_pair(1, COLOR_RED, COLOR_BLUE)

    t = "Hello World"
    h, w = stdscr.getmaxyx()
    x = w//2 - len(t)//2
    y = h//2

    textpad.rectangle(stdscr, 3, 3, h-3, w-3)

    stdscr.attron(color_pair(1))
    stdscr.addstr(y, x, t)
    stdscr.attroff(color_pair(1))

    stdscr.refresh()
    stdscr.getch()

wrapper(chrxp)
