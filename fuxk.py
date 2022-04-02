from math import *

def solve(func, aimv=0, es=0.01):
    '''简单的二分法求解。函数应为一元单调，定义在正值区间内。\n'''
    # 确定单调性，决定找根方向
    if (func(100) >= func(1)) ^ (func(10) < aimv):
        # 递增且在 10 处大于目标值，或者递减且在 10 处小于目标值，根需要在 10 之前找
        dire = 0.1
        v2 = 10
    else:
        # 递增且在 10 处大于目标值，或者递减且在 10 处小于目标值，根需要在 10 之后找
        dire = 10
        v2 = 10
    # 找初始区间
    count = 0
    while True:
        count += 1
        if count>=1024:
            raise Exception('貌似没根。')
        v1 = v2
        v2 = v1*dire
        if (func(v1)>aimv) ^ (func(v2)>aimv):
            break
        if abs(func(v1)-aimv)<es:
            return round(v1, 1-int(log(es,10)))
        if abs(func(v2)-aimv)<es:
            return round(v2, 1-int(log(es,10)))
    # 二分法
    fv1 = func(v1)
    fv2 = func(v2)
    count = 0
    while count < 1024:
        count += 1
        v0 = (v1+v2)/2
        fv0 = func(v0)
        if abs(fv0-aimv)<es:
            return round(v0, 1-int(log(es,10)))
        if (fv0>aimv) ^ (fv1>aimv):
            if abs(fv0-fv1) < es:
                # 已经小于指定精度，直接返回
                return round((v1+v0)/2, 1-int(log(es,10)))
            v2 = v0
            fv2 = fv0
        else:
            if abs(fv0-fv2) < es:
                return round((v2+v0)/2, 1-int(log(es,10)))
            v1 = v0
            fv1 = fv0
    raise Exception('二分不出来，呜呜。')

class G:
    '''简单的系统罢了。就是这样的一个形式，再高级的不想搞，至少作业用不到，诶嘿~：
             K*(t1*s + 1)(t2*s + 1)...(tn*s + 1)
    G(s) = ---------------------------------------
            s^v*(T1*s + 1)(T2*s + 1)...(Tm*s + 1)
创建时需要的参数为：
- K: 开环增益，为正数
- tlist: 分子tn列表，记得套中括号，记得化为+1系统
- Tlist: 分母Tm列表，记得套中括号，记得化为+1系统
- v: 系统型别，也就是分母s因子的次数，为整数，默认为1，且仅为1时，系统才能计算指标
- es: 计算精度，二分法求解时有用，默认0.001，即保留2位小数

这些成员为系统即时的特性函数：
- amp(o) : 精确的开环幅频特性
- amp_log(o) : 大致的开环幅频特性，使用对数dB形式表示（可近似认为是 20*lg(amp(o))）
- amp_fb : 精确的闭环幅频特性
- phase(o) : 精确的相频特性，角度值

这些成员为系统即时的特性值：
- Wc : 精确的剪切频率
- Wc_log : 使用大致的对数幅频特性得到的剪切频率
- gamma : 精确的相角裕度
- gamma_log : 使用大致剪切频率得到的相角裕度
- Wg : 穿越频率
- Kg : 精确的幅值裕度
- Kg_log : 对数表示的、使用大致对数幅频特性得到的幅值裕度，使用对数形式表示
- Mr_log : 使用gamma_log计算得到的相对谐振峰值
- sigmap_approx : 使用经验公式，用gamma_log计算得到的超调量
- ts_approx : 使用经验公式，用gamma_log和Wc_log计算得到的5%调节时间
你可以使用 showinfo 方法来查看这些值，
此外，你也可以用 showlogexp 方法输出幅频对数算式，拿一点过程分嘛

你甚至可以玩频域校正！
- correct1 : 相角裕度优先的超前校正
- correct2 : 剪切频率优先的超前校正
- correct3 : 相角裕度优先的迟后校正
- correct4 : 剪切频率优先的迟后校正
不确定如何使用的话，可以用 print(G.correct1.__doc__) 看参数怎么给。

由于运算过程中的精度统一问题，所以一定要验算！
此外，如果有需要，可以 update 刷新特性函数和特性值。

加油算他！
'''
    def __init__(self, K, tlist, Tlist, v=1, es=0.001):
        self.tau = tlist.copy()
        self.time = Tlist.copy()
        self.v = v
        self.k = K

        self.es = es

        self.amp = None
        self.amp_log = None
        self.amp_fb = None
        self.phase = None
        self.Wc = 0.0
        self.gamma = 0.0
        self.Wc_log = 0.0
        self.gamma_log = 0.0
        self.Wg = 0.0
        self.Kg = 0.0
        self.Kg_log = 0.0
        self.Wb = 0.0
        self.Mr_log = 0.0
        self.sigmap_approx = 0.0
        self.ts_approx = 0.0

        self.ready = False

        self.update()

    def __repr__(self):
        et = 1-int(log(self.es,10))
        def rnd(n):
            return round(n, et)

        num = ''
        den = ''

        if not self.tau:
            num = str(self.k)
        elif self.k==1:
            num = ''.join([f'({rnd(t)}*s + 1)' for t in self.tau])
        else:
            num = str(self.k) + '*' + ''.join([f'({rnd(t)}*s + 1)' for t in self.tau])

        if not self.time:
            if self.v==0:
                den = '1'
            elif self.v==1:
                den = 's'
            else:
                den = f's^{self.v}'
        else:
            if self.v==0:
                den = ''.join([f'({rnd(t)}*s + 1)' for t in self.time])
            elif self.v==1:
                den = 's*' + ''.join([f'({rnd(t)}*s + 1)' for t in self.time])
            else:
                den = f's^{self.v}*' + ''.join([f'({rnd(t)}*s + 1)' for t in self.time])

        dn1 = len(num)
        dn2 = len(den)
        slash = '-'*max(dn1,dn2)
        if dn1>dn2:
            den = ' '*((dn1-dn2)//2) + den
        else:
            num = ' '*((dn2-dn1)//2) + num
        return num+'\n'+slash+'\n'+den

    def __mul__(self, gg):
        if type(gg) is int or type(gg) is float:
            return G(self.k*gg, self.tau, self.time, self.v, self.es)
        elif type(gg) is G:
            return G(self.k*gg.k, self.tau+gg.tau, self.time+gg.time, self.v+gg.v, max(self.es,gg.es))

    def update(self):
        self.ready = False
        self.tau.sort(reverse = True)
        self.time.sort(reverse = True)

        if any([t<=0 for t in self.tau]) or any([t<0 for t in self.time]) or self.v!=1 or self.k<=0:
            print('返回的系统有点高级哦，没法自动算指标嗷:(')
            return

        amp_exp = 'lambda o:' + str(self.k) + '/o**' + str(self.v)
        amp_log_exp = 'lambda o:20*log(' + str(self.k) + '/o**' + str(self.v)
        phase_exp = 'lambda o:' + str(-90*self.v) + '+(0'
        fbamp_num_exp = '1'

        for t in self.tau:
            t1 = 1/t
            amp_exp += f'*sqrt(1+({str(t)}*o)**2)'
            amp_log_exp += f'*(o*{str(t)} if o>={str(t1)} else 1)'
            fbamp_num_exp += f'*(1+(o*{str(t)}*1j))'
            phase_exp += f'+atan(o*{str(t)})'

        fbamp_den_exp = fbamp_num_exp+'+o*1j'

        for t in self.time:
            t1 = 1/t
            amp_exp += f'/sqrt(1+({str(t)}*o)**2)'
            amp_log_exp += f'/(o*{str(t)} if o>={str(t1)} else 1)'
            fbamp_den_exp += f'*(1+(o*{str(t)}*1j))'
            phase_exp += f'-atan(o*{str(t)})'

        amp_log_exp += ',10)'
        phase_exp += ')*180/pi'
        amp_fb_exp = 'lambda o:abs('+fbamp_num_exp+')/abs('+fbamp_den_exp+')'

        self.amp = eval(amp_exp)
        self.amp_log = eval(amp_log_exp)
        self.amp_fb = eval(amp_fb_exp)
        self.phase = eval(phase_exp)

        try:
            self.Wc = solve(self.amp, 1, self.es)
            self.Wc_log = solve(self.amp_log, 0, self.es)
            self.Wg = solve(self.phase, -180, self.es)
        except:
            print('寄啦！看来数有点刁钻，算不出来嗷。')
            print('下面是报错信息：')
            raise

        try:
            self.Wb = solve(self.amp_fb, 0.5**0.5, self.es)
        except:
            self.Wb = -1

        self.gamma = 180+self.phase(self.Wc)
        self.gamma_log = 180+self.phase(self.Wc_log)
        self.Kg = 1/self.amp(self.Wg)
        self.Kg_log = -self.amp_log(self.Wg)

        self.Mr_log = 1/sin(self.gamma_log*pi/180)
        singamma1_1 = self.Mr_log-1
        if 34<self.gamma_log<90:
            self.sigmap_approx = 0.16 + 0.4*(singamma1_1)
        else:
            self.sigmap_approx = -1
        self.ts_approx = pi/self.Wc_log*(2 + 1.5*singamma1_1 + 2.5*singamma1_1**2)

        self.ready = True

    def correct1(self, phim):
        '''简单而粗糙的相角裕度优先的超前校正而已。使用对数，精度不保证。\n参数就是想要超前的相角角度值。\n'''
        if phim<0 or phim>80:
            print('phim有些刁难哦，去死一死吧。')
            return None
        if not self.ready:
            print('原系统我都没算出来，校正不了啊呜')
            return None
        sinphi = sin(phim*pi/180)
        alpha = (1+sinphi)/(1-sinphi)
        om = solve(self.amp_log, -10*log(alpha, 10), self.es)
        T = 1/om/sqrt(alpha)
        print('返回了一个校正系统嗷，你可以将它与原系统相乘~~')
        return G(1, [alpha*T], [T], 0)

    def correct2(self, Wc):
        '''简单而粗糙的剪切频率优先的超前校正，使用对数，精度不保证。\n参数就是预期的剪切频率。\n'''
        if not self.ready:
            print('原系统我都没算出来，校正不了啊呜')
            return None
        if Wc<self.Wc_log:
            print('Wc这么小，让我超前个鬼哦！')
            return None
        taramp = self.amp_log(Wc)
        alpha = 10**(taramp/-10)
        T = 1/Wc/sqrt(alpha)
        print('返回了一个校正系统嗷，你可以将它与原系统相乘~~')
        return G(1, [alpha*T], [T], 0)

    def correct3(self, gammaplus, tk=7):
        '''简单而粗糙的相角裕度优先的迟后校正，使用对数，精度不保证。\n参数1就是考虑到迟后影响的预期相角裕度角度值，\n参数2为依据迟后影响决定的tau系数，一般为5~10，默认为7。\n'''
        if not self.ready:
            print('原系统我都没算出来，校正不了啊呜')
            return None
        if gammaplus < self.gamma_log:
            print('这比原系统相角裕度还小嗷，真就反向校正呗。')
            return None
        oc = solve(self.phase, gammaplus-180, self.es)
        taramp = self.amp_log(oc)
        beta = 10**(taramp/20)
        tau = tk/oc
        print('返回了一个校正系统嗷，你可以将它与原系统相乘~~')
        return G(1, [tau], [beta*tau], 0)

    def correct4(self, Wc, tk=7):
        '''简单而粗糙的剪切频率优先的迟后校正，使用对数，精度不保证。\n参数1就是考虑到迟后影响的预期剪切频率，\n参数2为依据迟后影响决定的tau系数，一般为5~10，默认为7。\n'''
        if not self.ready:
            print('原系统我都没算出来，校正不了啊呜')
            return None
        if Wc > self.Wc_log:
            print('这比原系统剪切频率还大嗷，确定是用迟后校正嘛。')
            return None
        taramp = self.amp_log(oc)
        beta = 10**(taramp/20)
        tau = tk/oc
        print('返回了一个校正系统嗷，你可以将它与原系统相乘~~')
        return G(1, [tau], [beta*tau], 0)

    def showinfo(self):
        print()
        print(self)
        print()
        if self.ready:
            print('开环指标精确值们：')
            print('\tWc\t=', self.Wc, 'rad/s')
            print('\tWg\t=', self.Wg, 'rad/s')
            print('\tGamma\t=', str(self.gamma)+'°')
            print('\tKg\t=', self.Kg)
            print('开环指标对数简化运算得到的值：')
            print('\tWc\t=', self.Wc_log, 'rad/s')
            print('\tWg\t=', self.Wg, 'rad/s')
            print('\tGamma\t=', str(self.gamma_log)+'°')
            print('\t20lg Kg\t=', self.Kg_log, 'dB')
            print('闭环频域近似指标：')
            print(f'\tWb\t= {self.Wb} rad/s' if self.Wb>0 else '\t不存在截止频率')
            print('\tMr\t=', self.Mr_log*100, r'%')
            print('闭环时域指标（经验公式、近似指标）：')
            print('\tts\t=', self.ts_approx, 's')
            print(f'\tsigma_p\t= {self.sigmap_approx*100} %' if self.sigmap_approx>0 else '\t相角裕度超出范围，超调量无法用经验公式计算')
            print()
        else:
            print('我太弱小了，没有力量，算不了指标，呜呜。。')
            print()

    def showlogexp(self):
        if self.ready:
            taus = self.tau.copy()
            times = self.time.copy()
            et = 1-int(log(self.es,10))
            def rnd(n):
                return round(n, et)

            print('20lg A(w) =')
            logexp = ' 20('
            if self.k!=1:
                logexp += f' lg {self.k} -lg w '  # 反正只用来算1型系统，直接写上了。
            else:
                logexp += '-lg w '

            t0 = 0
            while taus or times:
                if taus and times:
                    if taus[0] > times[0]:
                        t = taus.pop(0)
                        pm = '+'
                    else:
                        t = times.pop(0)
                        pm = '-'
                elif taus:
                    t = taus.pop(0)
                    pm = '+'
                elif times:
                    t = times.pop(0)
                    pm = '-'

                t1 = 1/t
                print(logexp+f') , {rnd(t0)} < w < {rnd(t1)}')
                logexp += f'{pm}lg {rnd(t)}w '
                t0 = t1
            print(logexp+f') , w > {rnd(t0)}')

        else:
            print('这个系统不支持这项功能哦、')

print('算TM的系统校正！')
