from curses import *
from curses import textpad

_istr = []

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
                
def _txtbox(stdscr, y, xl, wl = 20, xpndx = False):
    wl += xl+2
    istr = ''
    if not xpndx:
        textpad.rectangle(stdscr, y, xl, y+2, wl)
        stdscr.addstr(y+1, xl+1, '')
        i = xl+1
        while i < wl:
            k = stdscr.getch()
            if k == KEY_ENTER or k in [10, 13]:
                break
            elif k == KEY_BACKSPACE:
                if i>xl+1:
                    stdscr.addstr(y+1, i-1, ' ')
                    stdscr.addstr(y+1, xl+1, istr[:-1])
                    istr = istr[:-1]
                    i-=1
                    stdscr.refresh()
            else:
                if i<wl-1:
                    stdscr.addstr(y+1, i, str(chr(k)))
                    istr+=str(chr(k))
                    stdscr.refresh()
                    i+=1

    else:
        # wait for update...
        textpad.rectangle(stdscr, y, xl, y+2, wl)
        stdscr.addstr(y+1, xl+1, '')
        i = xl+1
        while True:
            if i == wl:
                textpad.rectangle(stdscr, y, xl, y+2, wl+i-xl)
                stdscr.addstr(y+1, i+1, '')
            k = stdscr.getch()
            if k == KEY_ENTER or k in [10, 13]:
                break
            elif k == KEY_BACKSPACE:
                if i>xl+1:
                    stdscr.addstr(y+1, i-1, ' ')
                    stdscr.addstr(y+1, xl+1, istr[:-1])
                    istr = istr[:-1]
                    i-=1
                    stdscr.refresh()
            else:
                if i<wl-1:
                    stdscr.addstr(y+1, i, str(chr(k)))
                    istr+=str(chr(k))
                    stdscr.refresh()
                    i+=1
            


    _istr.append(istr)

                
wrapper(chrxp)
