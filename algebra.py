__brkts__ = (('[', ']'),
             ('{', '}'),
             ('(', ')'))
def msum( ex: str):
    opts = []
    optr = []
    c = ''
    for i in ex:
        if i == '+' or i == '-':
            opts.append(c)
            optr.append(i)
            c = ''
        else: c+=i
    else: opts.append(c)
    dtor, ntor = '', ''
    for i in opts:
        dtor += i[i.find('/')+1:]
    for j, i in enumerate(opts):
        if j < len(optr):
            ntor += opts[j][:opts[j].index('/')] + div(dtor+'/'+i[i.index('/')+1:])+optr[j]
    else: ntor += opts[j][:opts[j].index('/')] + div(dtor+'/'+i[i.index('/')+1:])
    return cmn(ntor)+'/'+cmn(dtor)


def __cmn__(ex : str):
    if '+' in ex or '-' in ex:
        opts = []
        optr = []
        c = ''
        for i in ex:
            if i=='+' or i=='-':
                opts.append(c)
                optr.append(i)
                c = ''
            else:
                c += i
        else:
            opts.append(c)
        cmnvar = ''
        lopts = []
        for i in opts[0]:
            ct = 1
            for j in range(1, len(opts)):
                for k in opts[j]:
                    if i is k:
                        ct += 1
                        break
            if ct is len(opts): cmnvar += i
            for i in opts: lopts.append(list(i))
            for i in range(len(lopts)):
                for j in cmnvar:
                    if j in lopts[i]:
                        lopts[i].remove(j)
                ab = ''
                for j in lopts[i]:
                    ab += j
                if ab == '': ab = '1'
                lopts[i] = ab
        rval = cmnvar+'('
        for i in range(len(optr)): rval += lopts[i]+optr[i]
        else: rval += lopts[len(lopts)-1] + ')'
        return rval
    else:
        return ex

def __div__(itm: str):
    if '/' in itm:
        l1, c = list(__cmn__(itm.split('/')[0])), list(__cmn__(itm.split('/')[0]))
        l2 = list(itm.split('/')[1])
        if '(' in l1 and ')' in l1: del l1[l1.index('('):l1.index(')')+1]
        if '(' in l2 and ')' in l2: del l2[l2.index('('):l2.index(')')+1]
        for i in c:
            if i in l2:
                l1.remove(i)
                l2.remove(i)
        a, b = '', ''
        if l1 < l2:
            for i in range(len(l2)):
                if i < len(l1):
                    a += l1[i]
                b += l2[i]
        else:
            for i in range(len(l1)):
                if i < len(l2):
                    b += l2[i]
                a += l1[i]
        if a is "": a = '1'
        if b is '': return a
        return a + '/' + b
    else:
        return itm


#For testing:
#print(__div__('has/as'))
