import tkinter as tk
from tkinter import *
from tkinter.ttk import *
from tkinter.filedialog import askopenfilename
from tkinter import ttk
import numpy as np
from tkinter import messagebox
from time import time, sleep
from random import choice, uniform
from math import sin, cos
from sys import modules
import random
import subprocess
import os
from tkinter import simpledialog

class Interface_livrable_c():
    def __init__(self):
        self.xx = 0
        self.yy = 0
        self.compte = 5
        self.Interface_livrable_c = tk.Tk()
        self.Interface_livrable_c.title("Démonstration Livrable C")
        self.Interface_livrable_c.geometry('250x250')
        self.Interface_livrable_c.configure(background='black')
        #button = tk.Button(self.Interface_livrable_c, text = 'Exécuter', command=self.compute_c,bg='brown',fg='white',width=17)
        #button.grid(row=1, column=0)
        self.compute_c()

    def compute_c(self):
        if self.compte > 0:
            if self.compte != 5:
                self.Interface_livrable_c.after(200, self.dec.destroy())
            self.dec = Label(self.Interface_livrable_c, text=self.compte,width=20,font=("bold", 20), background="black",foreground="green" , borderwidth=0)
            self.dec.place(x=120, y=100)
            self.compte -= 1
            self.Interface_livrable_c.after(1000, self.compute_c)
        elif self.compte == 0:
            self.Interface_livrable_c.after(200, self.dec.destroy())
            self.dec = Label(self.Interface_livrable_c, text='Exécution!',width=20,font=("bold", 18), background="black",foreground="green")
            self.dec.place(x=60, y=100)
            self.compte -= 1
            self.Interface_livrable_c.after(2000, self.compute_c)
        elif self.compte > - 200:
            if self.compte == -1:
                self.dec.destroy()
            self.compte -= 1
            DATA = np.random.randint(2)
            tk.Label(self.Interface_livrable_c, text=DATA, background="black",foreground="green").place(x=self.xx,y=self.yy)
            self.xx += 20
            if self.xx > 250:
                self.xx = 0
                self.yy += 20
            self.Interface_livrable_c.after(20, self.compute_c)
        else :
            tmp=subprocess.call("./noice_c")
            print("\nExecution completed... \n")
            self.Interface_livrable_c.after(100, self.quit)



    def quit(self):
       self.Interface_livrable_c.destroy()


class Interface_view():
    def __init__(self):
        self.Interface_view = tk.Tk()
        self.Interface_view.title("Vérification")
        self.Interface_view.geometry('520x470')

    def compute(self, txtFile):
        self.read(txtFile)

    def read(self, txtFile):
        f= open(txtFile,"r")
        if f.mode == 'r':
            contents =f.read()
            print(contents)
        self.T = tk.Text(self.Interface_view, height=26, width=63)
        self.T.pack()
        self.T.insert(tk.END, contents)

    def quit(self):
       self.Interface_view.destroy()


class Interface():

    def __init__(self):
       self.interface = tk.Tk()
       self.interface.title("AER8875 - Projet intégrateur IV [Groupe D]")
       self.interface.geometry('600x750')
       tab_control = ttk.Notebook(self.interface)
       self.tab1 = ttk.Frame(tab_control)
       self.tab2 = ttk.Frame(tab_control)
       self.tab30 = ttk.Frame(tab_control)
       tab_control.add(self.tab1, text='Entrees')
       tab_control.add(self.tab2, text='Resolution')
       tab_control.add(self.tab30, text='Sorties')
       tab_control.pack(expand=1, fill='both')
       ######################## menu bar ########################
       menubar = Menu(self.interface)
       filemenu = Menu(menubar, tearoff=0)
       filemenu.add_command(label="Enregistrer", command=self.saveNext)
       filemenu.add_command(label="Recommencer", command=self.re_start)
       filemenu.add_command(label="Aide", command=self.help)
       filemenu.add_separator()
       filemenu.add_command(label="Fermer", command=self.quit)
       menubar.add_cascade(label="Options", menu=filemenu)
       checking = Menu(menubar, tearoff=0)
       checking.add_command(label="Vérification de la paramétrisation", command=self.saveNext)
       menubar.add_cascade(label="Vérification", menu=checking)
       test_case = Menu(menubar, tearoff=0)
       #test_case.add_command(label="Livrable C", command=self.livrableC)
       test_case.add_command(label="Livrable C [option nd]", command=self.ND_errorMSG)
       menubar.add_cascade(label="Cas tests", menu=test_case)
       self.interface.config(menu=menubar)
        #############################################
       formTitle = Label(self.tab1, text="Interface utilisateur",width=20,font=("bold", 18))
       formTitle.place(x=185,y=20)
       self.check_1 = BooleanVar()
       self.check_2 = BooleanVar()
       self.check_3 = BooleanVar()
       self.check_1.set(False)
       self.check_2.set(False)
       self.check_3.set(False)
       self.Limiteur = StringVar()
       self.optionLISTE = StringVar()
       self.delimiter4 = 0
       self.delimiter5 = 0
       self.boolCP_chord = BooleanVar()
       self.boolResiduConv = BooleanVar()
       self.boolCLConv = BooleanVar()
       self.boolCDConv = BooleanVar()
       self.boolCL_Alpha = BooleanVar()
       self.boolCD_Alpha = BooleanVar()
       self.boolCM_Alpha = BooleanVar()
       self.boolCP_chord.set(False)
       self.boolResiduConv.set(False)
       self.boolCLConv.set(False)
       self.boolCDConv.set(False)
       self.boolCL_Alpha.set(False)
       self.boolCD_Alpha.set(False)
       self.boolCM_Alpha.set(False)
       self.chgt = 0
       self.window = False
       self.interface_2 = False
       self.partition = 0
       self.save = 0
    ######################## Variables for text file ########################
       self.SimName=StringVar()
       self.choiceMESH = StringVar()
       self.choiceMESH.set("none")
       self.su2FilePath=StringVar()
       self.NbrPartition=IntVar()
       self.alpha=0.0
       self.mach=0.0
       self.gamma=1.4
       self.TemporelMethod = StringVar()
       self.stgaeNbr = IntVar()
       self.cflVAR_str=StringVar()
       self.ResidualSmooth = StringVar()
       self.SpatialMethod = StringVar()
       self.SpatialMethodORDER = StringVar()
       self.grad = StringVar()
       self.iterMAX = DoubleVar()
       self.convergenceCriterea = DoubleVar()
       self.graph = []
       self.alpha_0=0.0
       self.alpha_inf=0.0

    ######################## Name the project ########################
       lbl = tk.Label(self.tab1, text="Nommer la simulation :   ")
       #lbl.grid(column=0,row=1)
       lbl.place(x=10,y=90)
       self.txt = tk.Entry(self.tab1, width=20)
       self.txt.insert(0, "CFDsimPI4")
       #txt.grid(column=1,row=1)
       self.txt.place(x=170,y=90)
    ######################## Delimiters ########################
       #delimiter1 = Label(self.tab1, text="------------------------------------------------------------------------------------------------------------------")
       #delimiter1.place(x=10,y=160)
       #delimiter2 = Label(self.tab1, text="------------------------------------------------------------------------------------------------------------------")
       #delimiter2.place(x=10,y=250)
    ######################## Choice for mesh ########################
       self.btn_meshChoice_1 = tk.Button(self.tab1, text="Importer un maillage",height=2, width=54, command=self.MeshChoiceIMPORT)
       self.btn_meshChoice_1.place(x=70,y=130)
    ######################## Option : Patitionnement ########################
       lbl3 = tk.Label(self.tab1, text="Partitionnement du maillage :   ")
       lbl3.place(x=10,y=260)
       var = StringVar()
       var.set("4")
       self.spin_partition = Spinbox(self.tab1, bd=3, wrap=True, from_=1, to=15, width=5,textvariable=var)
       self.spin_partition.place(x=210,y=260)
    ######################## Option : Patitionnement ########################
       lbl333 = tk.Label(self.tab1, text="Nombre de coeur :   ")
       lbl333.place(x=10,y=300)
       var = StringVar()
       var.set("4")
       self.spin_partition1 = Spinbox(self.tab1, bd=3, wrap=True, from_=1, to=15, width=5,textvariable=var)
       self.spin_partition1.place(x=150,y=300)
    ######################## Option : alpha ########################
       lbl2 = tk.Label(self.tab1, text="Angle d'attaque :")
       lbl2.place(x=10,y=310+40)
       lbl22 = tk.Label(self.tab1, text="(\u03B1 [deg])")
       lbl22.place(x=185,y=310+40)
       self.entre1 = tk.Entry(self.tab1, width=5)
       self.entre1.insert(0, 0.0)
       self.entre1.place(x=135,y=310+40)
    ######################## Option : Mach Nbr ########################
       lbl2 = tk.Label(self.tab1, text="Nombre de Mach :")
       lbl2.place(x=10,y=360+40)
       self.entre2 = tk.Entry(self.tab1, width=5)
       self.entre2.insert(0, 1.0)
       self.entre2.place(x=135,y=360+40)
    ######################## Option : Gamma ########################
       lbl2 = tk.Label(self.tab1, text="Gamma :")
       lbl2.place(x=10,y=410+40)
       self.entre3 = tk.Entry(self.tab1, width=5)
       self.entre3.insert(0, self.gamma)
       self.entre3.place(x=80,y=410+40)
    ######################## LOGO ########################
       image = tk.PhotoImage(file="logo.png")
       image = image.subsample(6, 6)
       label = tk.Label(self.tab1,image=image)
       label.place(x=300,y=310)
################################################################################################################################################################
       self.spin_stage = Spinbox(self.tab2, bd=3, wrap=True, from_=1, to=5, width=5)
       formTitle_2 = Label(self.tab2, text="Integration spatiale",width=20,font=("bold", 14))
       formTitle_2.place(x=10,y=30)
       self.formTitle_3 = Label(self.tab2, text="Integration temporelle",width=20,font=("bold", 14))
       self.formTitle_3.place(x=10,y=210)
    ######################## Option : Schema spatial ########################
       OptionList1 = [
       "Roe",
       "AUSM"]
       self.SpatialMethod.set("Choisir une méthode d'intégration spatiale")
       self.opt11 = tk.OptionMenu(self.tab2, self.SpatialMethod, *OptionList1)
       self.opt11.config(width=50, font=('Helvetica', 12))
       self.opt11.place(x=50,y=70)
       OptionList10 = [
       "Schema spatial d'ordre 1",
       "Schema spatial d'ordre 2"]
       self.SpatialMethodORDER.set("Ordre d'intégration")
       self.opt1 = tk.OptionMenu(self.tab2, self.SpatialMethodORDER, *OptionList10)
       self.opt1.config(width=20, font=('Helvetica', 12))
       self.opt1.place(x=50,y=110)
       self.SpatialMethodORDER.trace("w", self.SecOrdreEffect)
       self.SpatialMethod.trace("w", self.SpatialeBTN)
 ######################## Option : Schema temporel ########################
       OptionList2 = [
       "Euler-explicite",
       "Euler-implicite",
       "Runge-Kutta",
       "Hybride"]
       self.TemporelMethod.set("Choisir une méthode d'intégration temporelle")
       self.opt10 = tk.OptionMenu(self.tab2, self.TemporelMethod, *OptionList2)
       self.opt10.config(width=50, font=('Helvetica', 12))
       self.opt10.place(x=50,y=250)
       self.TemporelMethod.trace("w", self.TemporelMethod_BTN)
    ######################## Option : cfl ########################
       self.lbl1 = tk.Label(self.tab2, text="Coefficient de Courant–Friedrichs–Lewy :")
       self.lbl1.place(x=10,y=310)
       self.lbl11 = tk.Label(self.tab2, text="(CFL)")
       self.lbl11.place(x=330,y=310)
       self.txt1 = tk.Entry(self.tab2, width=5)
       self.txt1.insert(0, 1.0)
       self.txt1.place(x=285,y=310)
    ######################## Option : Smoothing ########################
       self.lbl111 = tk.Label(self.tab2, text="Lissage du résidu :")
       self.lbl111.place(x=10,y=360)
       OptionList1 = [
       "Yes",
       "No"]
       self.ResidualSmooth.set("Yes or No")
       self.opt2 = tk.OptionMenu(self.tab2, self.ResidualSmooth, *OptionList1)
       self.opt2.config(width=8, font=('Helvetica', 10))
       self.opt2.place(x=140,y=355)
       self.ResidualSmooth.trace("w", self.Lissage_BTN)
################################################################################################################################################################
       self.formTitle_4 = Label(self.tab2, text="Convergence",width=20,font=("bold", 14))
       self.formTitle_4.place(x=10,y=450)

       self.label3 = tk.Label(self.tab2, text="Nombre d'iteration maximale :")
       self.label3.place(x=10,y=490)
       self.text2 = tk.Entry(self.tab2, width=10)
       self.text2.insert(0, 100000.0)
       self.text2.place(x=240,y=490)

       self.label4 = tk.Label(self.tab2, text="Critere de convergence :")
       self.label4.place(x=10,y=530)
       self.text3 = tk.Entry(self.tab2, width=8)
       self.text3.insert(0, 1e-14)
       self.text3.place(x=200,y=530)

       btn1_info = tk.Button(self.tab2, text="?",bg='gold',fg='black',width=1,height=1)
       btn1_info.place(x=270,y=527)
       text_info = ('Notez que la convergence est verifiee '
       'en faisant la somme des carree des residues '
       'pondere sur les volumes ')
       self.waittime = 500     #miliseconds
       self.wraplength = 180   #pixels
       self.widget = btn1_info
       self.text_info = text_info
       self.widget.bind("<Enter>", self.enter)
       self.widget.bind("<Leave>", self.leave)
       self.widget.bind("<ButtonPress>", self.leave)
       self.id = None
       self.tw = None

######################## Option : Graph to generate ########################
       formTitle1 = Label(self.tab30, text="Indiquer les graphiques à générer",font=("bold", 15))
       formTitle1.place(x=125,y=53)
       lbl5 = tk.Label(self.tab30, text="Convergence du résidu :   ")
       lbl5.place(x=10,y=100)
       chk1 = Checkbutton(self.tab30, var=self.boolResiduConv)
       chk1.place(x=500,y=100)

       lbl8 = tk.Label(self.tab30, text="Distribution de pression sur le profil :   ")
       lbl8.place(x=10,y=130)
       chk2 = Checkbutton(self.tab30, var=self.boolCP_chord)
       chk2.place(x=500,y=130)

       lbl6 = tk.Label(self.tab30, text="Convergence du coefficient de portance :   ")
       lbl6.place(x=10,y=160)
       chk3 = Checkbutton(self.tab30, var=self.boolCLConv)
       chk3.place(x=500,y=160)

       lbl7 = tk.Label(self.tab30, text="Convergence du coefficient de traînée :   ")
       lbl7.place(x=10,y=190)
       chk4 = Checkbutton(self.tab30, var=self.boolCDConv)
       chk4.place(x=500,y=190)

       lbl8 = tk.Label(self.tab30, text="Coefficient de portance en fonction de l'angle d'attaque :   ")
       lbl8.place(x=10,y=220)
       chk5 = Checkbutton(self.tab30, var=self.boolCL_Alpha)
       chk5.place(x=500,y=220)
       self.boolCL_Alpha.trace("w", self.Cx_Alpha)

       lbl9 = tk.Label(self.tab30, text="Coefficient de traînée en fonction de l'angle d'attaque :   ")
       lbl9.place(x=10,y=250)
       chk6 = Checkbutton(self.tab30, var=self.boolCD_Alpha)
       chk6.place(x=500,y=250)
       self.boolCD_Alpha.trace("w", self.Cx_Alpha)

       lbl88 = tk.Label(self.tab30, text="Coefficient de moment de tengage en fonction de l'angle d'attaque :   ")
       lbl88.place(x=10,y=280)
       chk55 = Checkbutton(self.tab30, var=self.boolCM_Alpha)
       chk55.place(x=500,y=280)
       self.boolCM_Alpha.trace("w", self.Cx_Alpha)
    ######################## Position moment et corde et surface ########################
       formTitle_30 = Label(self.tab30, text="Précisions pour le calcul des coefficients aérodynamiques",font=("bold", 13))
       formTitle_30.place(x=10,y=360)

       self.Sref = tk.Entry(self.tab30, width=18)
       self.Sref.insert(0, 1.0)
       self.Sref.place(x=100,y=400)
       lbl300 = tk.Label(self.tab30, text="Surface de référence : ")
       lbl300.place(x=10,y=400)

       self.Cref = tk.Entry(self.tab30, width=18)
       self.Cref.insert(0, 1.0)
       self.Cref.place(x=100,y=430)
       lbl300 = tk.Label(self.tab30, text="Corde de référence : ")
       lbl300.place(x=10,y=430)

       self.xref = tk.Entry(self.tab30, width=4)
       self.xref.insert(0, 0.25)
       self.xref.place(x=310,y=460)
       self.yref = tk.Entry(self.tab30, width=4)
       self.yref.insert(0, 0.00)
       self.yref.place(x=350,y=460)
       self.zref = tk.Entry(self.tab30, width=4)
       self.zref.insert(0, 0.00)
       self.zref.place(x=390,y=460)
       lbl300 = tk.Label(self.tab30, text="Position pour le calcul de moment [x,y,z] : ")
       lbl300.place(x=10,y=460)
    ######################## BoutonsDuBas ########################
       button = tk.Button(self.tab1, text = 'Abandoner', command=self.quit,bg='brown',fg='white',width=17)
       button.place(x=60,y=600)
       button2 = tk.Button(self.tab1, text = 'Sauvegarder', command=self.saveNext,bg='brown',fg='white',width=17)
       button2.place(x=360,y=600)

       button = tk.Button(self.tab2, text = 'Abandoner', command=self.quit,bg='brown',fg='white',width=17)
       button.place(x=60,y=600)
       button2 = tk.Button(self.tab2, text = 'Sauvegarder', command=self.saveNext,bg='brown',fg='white',width=17)
       button2.place(x=360,y=600)

       button11 = tk.Button(self.tab30, text = 'Abandoner', command=self.quit,bg='brown',fg='white',width=17)
       button11.place(x=60,y=600)
       button22 = tk.Button(self.tab30, text = 'Sauvegarder', command=self.saveNext,bg='brown',fg='white',width=17)
       button22.place(x=360,y=600)

       style = Style()
       style.configure('TButton', font =('calibri', 20, 'bold'),borderwidth = '4')
       button3 = tk.Button(self.tab1, text = 'Lancer la simulation', command=self.computeCODE,bg='brown',fg='white',width=55)
       button3.place(x=60,y=650)
       button33 = tk.Button(self.tab30, text = 'Lancer la simulation', command=self.computeCODE,bg='brown',fg='white',width=55)
       button33.place(x=60,y=650)
       button333 = tk.Button(self.tab2, text = 'Lancer la simulation', command=self.computeCODE,bg='brown',fg='white',width=55)
       button333.place(x=60,y=650)

    #####################################################
    ######################## RUN ########################
    #####################################################
       self.interface.mainloop()
############################################################################################################################
#######################################################   METHODS   ########################################################
############################################################################################################################
    def seeTXT(self, *arg):
        if self.interface_2 == True:
            try:
                self.check.quit()
            except:
                self.interface_2 = True
        self.interface_2 = True
        self.window = True
        self.check = Interface_view()
        txtFILE = self.SimName + ".txt"
        self.check.compute(txtFILE)

    def ND_errorMSG(self):
        messagebox.showerror("Erreur", "Cette option n'est pas encore disponible ou a été désactiver")

    def re_start(self, *arg):
        #time.sleep(10)
        maker = Constructor()
        maker.restart(self)


    def help(self, *arg):
        print("Ouverture du guide utilisateur...")

    def enter(self, event=None):
        self.schedule()
    def leave(self, event=None):
        self.unschedule()
        self.hidetip()
    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(self.waittime, self.showtip)
    def unschedule(self):
        id = self.id
        self.id = None
        if id:
            self.widget.after_cancel(id)
    def showtip(self, event=None):
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        # creates a toplevel window
        self.tw = tk.Toplevel(self.widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tw, text=self.text_info, justify='left',
        background="#ffffff", relief='solid', borderwidth=1,
        wraplength = self.wraplength)
        label.pack(ipadx=1)
    def hidetip(self):
        tw = self.tw
        self.tw= None
        if tw:
            tw.destroy()
            # testing ...
            #
    def SpatialeBTN(self, *arg):
        if self.SpatialMethod.get()=="AUSM":
            self.opt11.destroy()
            OptionList1 = [
            "Roe",
            "AUSM"]
            self.SpatialMethod.set("Choisir une méthode d'intégration spatiale")
            self.opt11 = tk.OptionMenu(self.tab2, self.SpatialMethod, *OptionList1)
            self.opt11.config(width=50, font=('Helvetica', 12))
            self.opt11.place(x=50,y=70)
            messagebox.showinfo("Note", "Cette option n'est pas encore disponible ou a été désactiver")
        else:
            self.opt11.config(width=50, font=('Helvetica', 12),bg='LightBlue3')

    def MeshChoiceIMPORT(self, *arg):
        self.btn_meshChoice_1.destroy()
        self.btn_meshChoice_1 = tk.Button(self.tab1, text="Importer un maillage",height=2, width=54,command=self.MeshChoiceIMPORT,bg='LightBlue3')
        self.btn_meshChoice_1.place(x=70,y=130)
        self.choiceMESH = "IMPORT"
    ######################## Import SU2 file ########################
        try:
            self.lbl2.destroy()
            self.FileSearch_btn.destroy()
            self.lbl20.destroy()
            self.lbl2 = tk.Label(self.tab1, text="Choisir un maillage")
            self.lbl2.place(x=10,y=200)
            self.photo = PhotoImage(file = r"loupe.png")
            self.photoimage = self.photo.subsample(20, 20)
            self.FileSearch_btn = Button(self.tab1, image = self.photoimage, command = lambda:self.get_filePath())
            self.FileSearch_btn.place(x=145,y=190)
            self.lbl20 = tk.Label(self.tab1, text="...........................................su2",width=42)
            self.lbl20.place(x=230,y=200)
        except:
            self.lbl2 = tk.Label(self.tab1, text="Choisir un maillage")
            self.lbl2.place(x=10,y=200)
            self.photo = PhotoImage(file = r"loupe.png")
            self.photoimage = self.photo.subsample(20, 20)
            self.FileSearch_btn = Button(self.tab1, image = self.photoimage, command = lambda:self.get_filePath())
            self.FileSearch_btn.place(x=145,y=190)
            self.lbl20 = tk.Label(self.tab1, text="...........................................su2",width=42)
            self.lbl20.place(x=230,y=200)

    def ComputePartitions(self):
        self.btn_partition = tk.Button(self.tab1, text="Partitionner le maillage",height=2,command=self.ComputePartitions,bg='LightBlue3')
        self.btn_partition.place(x=310,y=250)
        calll = "./bin/Execute1 " + self.txt.get() + " " + self.su2FilePath + " " + self.spin_partition.get() + " GUI"
        os.system(calll)
        self.partition = 1
        print("\nPartitionning completed... \n")
        messagebox.showinfo(title="Opération réussi", message="Partitionnement terminé")

    def SecOrdreEffect(self, *arg):
        self.opt1.destroy()
        OptionList10 = [
        "Schema spatial d'ordre 1",
        "Schema spatial d'ordre 2"]
        self.opt1 = tk.OptionMenu(self.tab2, self.SpatialMethodORDER, *OptionList10)
        self.opt1.config(width=20, font=('Helvetica', 12))
        self.opt1.place(x=50,y=110)
        self.opt1.config(width=20, font=('Helvetica', 12),bg='LightBlue3')

    def TemporelMethod_BTN(self, *arg):
        self.opt10.config(width=50, font=('Helvetica', 12),bg='LightBlue3')
        if self.TemporelMethod.get()=="Runge-Kutta":
            self.opt10.config(width=20, font=('Helvetica', 12))
            self.opt10.place(x=50,y=250)
            self.stage = tk.Label(self.tab2, text="Nombre de stage :   ")
            self.stage.place(x=320,y=255)
            self.spin_stage = Spinbox(self.tab2, bd=3, wrap=True, from_=1, to=5, width=5)
            self.spin_stage.place(x=450,y=255)
            self.check_3.set(True)
        elif self.TemporelMethod.get()=="Hybride":
            self.opt10.config(width=20, font=('Helvetica', 12))
            self.opt10.place(x=50,y=250)
            self.stage = tk.Label(self.tab2, text="Nombre de stage :   ")
            self.stage.place(x=320,y=255)
            self.spin_stage = Spinbox(self.tab2, bd=3, wrap=True, from_=1, to=5, width=5)
            self.spin_stage.place(x=450,y=255)
            self.check_3.set(True)
        else:
            if self.check_3.get() == True:
                self.spin_stage.destroy()
                self.stage.destroy()
            self.opt10.config(width=50, font=('Helvetica', 12))
            self.opt10.place(x=50,y=250)
            self.check_3.set(False)


    def Lissage_BTN(self, *arg):
        self.opt2.config(width=8, font=('Helvetica', 10),bg='LightBlue3')

    def Cx_Alpha(self, *arg):
        if self.boolCL_Alpha.get() == True:
            self.boolx_Alpha = True
        elif self.boolCD_Alpha.get() == True:
            self.boolx_Alpha = True
        elif self.boolCM_Alpha.get() == True:
            self.boolx_Alpha = True
        else :
            self.boolx_Alpha = False

        if self.boolx_Alpha==True:
            if self.chgt==0:
                self.self.lbl10 = tk.Label(self.tab30, text="Définir la plage de valeur d'angle d'attaque pour l'analyse :")
                self.self.lbl10.place(x=50,y=300)

                self.self.lbl11 = tk.Label(self.tab30, text="De :")
                self.self.lbl11.place(x=70,y=320)
                self.self.lbl12 = tk.Label(self.tab30, text="[deg]")
                self.self.lbl12.place(x=160,y=320)
                self.txt10 = tk.Entry(self.tab30, width=5)
                self.txt10.insert(0, 0.0)
                self.txt10.place(x=100,y=320)

                self.self.lbl13 = tk.Label(self.tab30, text="À :")
                self.self.lbl13.place(x=70,y=340)
                self.self.lbl14 = tk.Label(self.tab30, text="[deg]")
                self.self.lbl14.place(x=160,y=340)
                self.txt11 = tk.Entry(self.tab30, width=5)
                self.txt11.insert(0, 10.0)
                self.txt11.place(x=100,y=340)

                self.chgt = 1
        else :
            if self.chgt == 1:
                self.self.lbl10.destroy()
                self.self.lbl11.destroy()
                self.self.lbl12.destroy()
                self.txt10.destroy()
                self.self.lbl13.destroy()
                self.self.lbl14.destroy()
                self.txt11.destroy()

                self.chgt = 0

    def get_filePath(self):
        filePath = askopenfilename(parent=self.tab1, title='Choose an SU2 file for simulation',filetypes =[('SU2 Files', '*.su2')])
        if filePath is not None:
            self.lbl20.destroy()
            self.su2FilePath = filePath
            pos = 0
            for char in filePath:
                pos += 1
                if char == "/":
                    posSave = pos
                if char == ".":
                    posSave2 = pos
            try:
                self.lbl20.destroy()
                self.filepath_appear.destroy()
                show = filePath[posSave : posSave2+3]
                self.lbl20 = tk.Label(self.tab1, text="...........................................su2",bg="SlateGray1",width=42)
                self.filepath_appear = tk.Label(self.tab1, text=show,bg="SlateGray1",width=42)
                self.filepath_appear.place(x=230,y=200)
                ######################## Partionner le maillage ########################
                self.btn_partition = tk.Button(self.tab1, text="Partitionner le maillage",height=2,command=self.ComputePartitions)
                self.btn_partition.place(x=310,y=250)
            except:
                show = filePath[posSave : posSave2+3]
                self.lbl20 = tk.Label(self.tab1, text="...........................................su2",bg="SlateGray1",width=42)
                self.filepath_appear = tk.Label(self.tab1, text=show,bg="SlateGray1",width=42)
                self.filepath_appear.place(x=230,y=200)
                ######################## Partionner le maillage ########################
                self.btn_partition = tk.Button(self.tab1, text="Partitionner le maillage",height=2,command=self.ComputePartitions)
                self.btn_partition.place(x=310,y=250)



    def on_spinbox_change_partition(self):
        self.NbrPartition = self.spin_partition.get()

    def on_spinbox_change_stage(self):
        self.stgaeNbr = self.spin_stage.get()


    def computeCODE(self, *arg):
        if self.save == 0 :
            messagebox.showerror("Erreur", 'Pour lancer la simulation, il faut minimalement'
                ' importer un maillage, le partitionner et sauvegarder le paramétrage')
        else :
            coeur = self.spin_partition1.get()
            commande = "mpirun -n " + coeur + " ./bin/Execute2"
            ExecCommand = simpledialog.askstring(title=" ", prompt="Confirmer la commande d'exécution:", initialvalue=commande)
            os.system(ExecCommand)
            messagebox.showinfo(title="Opération réussi", message="Exécution de la simulation terminée")

    def saveNext(self):
        if self.partition == 0 :
            messagebox.showerror("Erreur", "Pour sauvegarder le paramétrage et générer le fichier de paramètres,"
                ' il faut minimalement importer un maillage et le partitionner')

        else :
            self.alpha = self.entre1.get()
            self.mach = self.entre2.get()
            self.gamma = self.entre3.get()
            self.cflVAR_str = self.txt1.get()
            self.convergenceCriterea = self.text3.get()
            self.iterMAX = self.text2.get()
            self.on_spinbox_change_partition()
            self.on_spinbox_change_stage()
            self.SimName = self.txt.get()
            writer = InputFILE(self)
            writer.compute(self)
            self.seeTXT()
            self.save = 1

    def quit(self):
        try:
            if self.window == True :
                self.check.quit()
            self.interface.destroy()
        except:
            self.interface.destroy()

    def livrableC(self):
        self.quit()
        self.livrableC_effect()

    def livrableC_effect(self, *arg):
        self.livrableC = Interface_livrable_c()

# =========================================================================================================== #
class Constructor():
    def compute(self):
        app = Interface()
    def restart(self, interfac_e):
        interfac_e.quit()
        self.compute()

class InputFILE():
    def __init__(self, Interface):
        self.fileNAME = Interface.SimName + ".txt"


    def compute(self, Interface):
        ordreSpatial = Interface.SpatialMethodORDER.get()
        if ordreSpatial=="Schema spatial d'ordre 2":
            ordreSpatial = 2
        else:
            ordreSpatial = 1
        graphss = [Interface.boolResiduConv.get(), Interface.boolCP_chord.get(), Interface.boolCLConv.get(), Interface.boolCDConv.get(), Interface.boolCL_Alpha.get(), Interface.boolCD_Alpha.get(), Interface.boolCM_Alpha.get()]
        try:
            f = open(self.fileNAME)
            # Do something with the file
            answer = messagebox.askyesno("Déjà une paramètrisation à ce nom","Écraser la dernière paramétrisation?")
            if answer == True:
                print('Écraser le fichier actuel')
                txtFILE = open(self.fileNAME,"w")
                txtFILE.write("------------------------------------------------------------\n "
                                "                   Input File                             \n"
                                "------------------------------------------------------------\n")
                txtFILE.write("\n")
                txtFILE.write(str("SimName : "))
                txtFILE.write(str(Interface.SimName))
                txtFILE.write("\n")
                txtFILE.write(str("choiceMESH : "))
                txtFILE.write(str(Interface.choiceMESH))
                txtFILE.write("\n")
                txtFILE.write(str("su2FilePath : "))
                txtFILE.write(str(Interface.su2FilePath))
                txtFILE.write("\n")
                txtFILE.write(str("NbrPartition : "))
                txtFILE.write(str(Interface.NbrPartition))
                txtFILE.write("\n")
                txtFILE.write(str("alpha : "))
                txtFILE.write(str(Interface.alpha))
                txtFILE.write("\n")
                txtFILE.write(str("mach : "))
                txtFILE.write(str(Interface.mach))
                txtFILE.write("\n")
                txtFILE.write(str("gamma : "))
                txtFILE.write(str(Interface.gamma))
                txtFILE.write("\n")
                txtFILE.write(str("TemporelMethod : "))
                txtFILE.write(str(Interface.TemporelMethod.get()))
                txtFILE.write("\n")
                txtFILE.write(str("stgaeNbr : "))
                txtFILE.write(str(Interface.stgaeNbr))
                txtFILE.write("\n")
                txtFILE.write(str("cflVAR_str : "))
                txtFILE.write(str(Interface.cflVAR_str))
                txtFILE.write("\n")
                txtFILE.write(str("ResidualSmooth : "))
                txtFILE.write(str(Interface.ResidualSmooth.get()))
                txtFILE.write("\n")
                txtFILE.write(str("SpatialMethod : "))
                txtFILE.write(str(Interface.SpatialMethod.get()))
                txtFILE.write("\n")
                txtFILE.write(str("SpatialMethodORDER : "))
                txtFILE.write(str(ordreSpatial))
                txtFILE.write("\n")
                txtFILE.write(str("grad : "))
                txtFILE.write(str(Interface.grad.get()))
                txtFILE.write("\n")
                txtFILE.write(str("iterMAX : "))
                txtFILE.write(str(Interface.iterMAX))
                txtFILE.write("\n")
                txtFILE.write(str("convergenceCriterea : "))
                txtFILE.write(str(Interface.convergenceCriterea))
                txtFILE.write("\n")
                txtFILE.write(str("graph : "))
                txtFILE.write(str(graphss))
                txtFILE.write("\n")
                txtFILE.write(str("alpha_0 : "))
                txtFILE.write(str(Interface.alpha_0))
                txtFILE.write("\n")
                txtFILE.write(str("alpha_inf : "))
                txtFILE.write(str(Interface.alpha_inf))
                txtFILE.write("\n")
                txtFILE.close()
            else:
                print('Créer un nouveau fichier')
                self.fileNAME = self.fileNAME + '_1'
                print(self.fileNAME)
                self.compute(Interface)
        except IOError:
            txtFILE = open(self.fileNAME,"w")
            txtFILE.write("------------------------------------------------------------\n "
                            "                   Input File                             \n"
                            "------------------------------------------------------------\n")
            txtFILE.write("\n")
            txtFILE.write(str("SimName : "))
            txtFILE.write(str(Interface.SimName))
            txtFILE.write("\n")
            txtFILE.write(str("choiceMESH : "))
            txtFILE.write(str(Interface.choiceMESH))
            txtFILE.write("\n")
            txtFILE.write(str("su2FilePath : "))
            txtFILE.write(str(Interface.su2FilePath))
            txtFILE.write("\n")
            txtFILE.write(str("NbrPartition : "))
            txtFILE.write(str(Interface.NbrPartition))
            txtFILE.write("\n")
            txtFILE.write(str("alpha : "))
            txtFILE.write(str(Interface.alpha))
            txtFILE.write("\n")
            txtFILE.write(str("mach : "))
            txtFILE.write(str(Interface.mach))
            txtFILE.write("\n")
            txtFILE.write(str("gamma : "))
            txtFILE.write(str(Interface.gamma))
            txtFILE.write("\n")
            txtFILE.write(str("TemporelMethod : "))
            txtFILE.write(str(Interface.TemporelMethod.get()))
            txtFILE.write("\n")
            txtFILE.write(str("stgaeNbr : "))
            txtFILE.write(str(Interface.stgaeNbr))
            txtFILE.write("\n")
            txtFILE.write(str("cflVAR_str : "))
            txtFILE.write(str(Interface.cflVAR_str))
            txtFILE.write("\n")
            txtFILE.write(str("ResidualSmooth : "))
            txtFILE.write(str(Interface.ResidualSmooth.get()))
            txtFILE.write("\n")
            txtFILE.write(str("SpatialMethod : "))
            txtFILE.write(str(Interface.SpatialMethod.get()))
            txtFILE.write("\n")
            txtFILE.write(str("SpatialMethodORDER : "))
            txtFILE.write(str(ordreSpatial))
            txtFILE.write("\n")
            txtFILE.write(str("grad : "))
            txtFILE.write(str(Interface.grad.get()))
            txtFILE.write("\n")
            txtFILE.write(str("iterMAX : "))
            txtFILE.write(str(Interface.iterMAX))
            txtFILE.write("\n")
            txtFILE.write(str("convergenceCriterea : "))
            txtFILE.write(str(Interface.convergenceCriterea))
            txtFILE.write("\n")
            txtFILE.write(str("graph : "))
            txtFILE.write(str(graphss))
            txtFILE.write("\n")
            txtFILE.write(str("alpha_0 : "))
            txtFILE.write(str(Interface.alpha_0))
            txtFILE.write("\n")
            txtFILE.write(str("alpha_inf : "))
            txtFILE.write(str(Interface.alpha_inf))
            txtFILE.write("\n")
            txtFILE.close()
        finally:
            print('Ecriture terminee...')



maker = Constructor()
maker.compute()
