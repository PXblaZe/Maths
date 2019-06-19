from curses import *
from curses import textpad
def chrxp(stdscr):
    curs_set(0)
    init_pair(1, COLOR_RED, COLOR_BLUE)
    init_pair(2, COLOR_BLACK, COLOR_WHITE)
    mousemask(1)
    while 1:
        stdscr.clear()
        t = "Hello World"
        h, w = stdscr.getmaxyx()
        x = w//2 - len(t)//2
        y = h//2
        k = h-4
        l = (w-3)//2

        textpad.rectangle(stdscr, 3, 3, h-3, w-3)
        stdscr.attron(color_pair(2))
        stdscr.addstr(k, l, "Exit")
        stdscr.attroff(color_pair(2))


        stdscr.attron(color_pair(1))
        stdscr.addstr(y, x, t)
        stdscr.attroff(color_pair(1))

        stdscr.refresh()

        ky = stdscr.getch()
        if ky == KEY_MOUSE:
            _, a, b, _, _ = getmouse()
            stdscr.addstr(0, 0,str(a)+' '+str(b))
            stdscr.refresh()
            if a>l-1 and a<(l+4) and b == k:
                break
wrapper(chrxp)
