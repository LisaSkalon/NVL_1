import numpy as np
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
import argparse

def parse_inputs():
     parser = argparse.ArgumentParser(description='Parser for GMS-track builder')
     parser.add_argument('-i', '--input' , help='Input file' , metavar='Str',
                    type=str, required=True)
     
     parser.add_argument('-s', '--string', help = 'String, printed before matrix', metavar = 'Str', type = str, default = 'Final matrix:')
     parser.add_argument('-g', '--gap', help = 'Score for gap', metavar = 'Int', type = int, default = -5)
     parser.add_argument('-p', '--print_flag', help='Print seq1 and seq 2 for me!', action='store_true')
     parser.add_argument('-n', '--name', help = 'File name', metavar = 'Str', type = str, default = 'our_output')
     args = parser.parse_args()
     return args.input, args.string, args.gap, args.print_flag, args.name


def levenstein(gapp, print_flag, seq1, seq2):
    # На всякий случай делaем регистр большим (в блосуме используются именно
    # заглавные буквы).
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    n = len(seq1)
    m = len(seq2)
    # Заполняем матрицу для поиска кратчайшего пути нулями.
    Wmin = np.zeros((n+1, m+1))
    # Заполняем нулями матрицы для сдвигов по вертикали и по горизонтали.
    FromV = np.zeros((n+1, m+1))
    FromG = np.zeros((n+1, m+1))
    # Создаем пустые списки, куда будем записывать выравненные 
    # последовательности.
    s1=[]
    s2=[]
    a = 0
    b = 0
    # Создаем строки, куда запишем выравненные последовательности.
    al_seq1 = ''
    al_seq2 = ''
    # Вводим значение штрафа за гэп
    gap = gapp

    # Проходим циклом по ячейкам матрицы. При этом горизонтальные и 
    # вертикальные переходы дают штраф в -5 (поскольку означают гэпы), а 
    # диагональные - соответствующий штраф из матрицы blosum62. В каждый
    # проход цикла заполняем соответствующую ячейку матрицы
    # значением штрафа наиболее выгодного варианта(самый большой - самый
    # выгодный). Нулевую горизонталь и
    # вертикаль заполняем вручную цифрами от 0 до n или m - поскольку мы 
    # точно знаем ценность вертикальных/горизонтальных переходов(это ценность
    # гэпов).
 
        
    for i in range(1,n+1):
        for j in range(1,m+1):
            
                Wmin[i, 0] = (-5)*i
                Wmin[0, j] = (-5)*j
                
                a1 = Wmin[i, j-1] + gap
                a2 = Wmin[i-1, j] + gap
                # В матрице blosum хранятся только значения штрафов для (x,y).
                # С помощью переворачивания получаем нужное значение для (y,x).
                pair = seq1[i-1], seq2[j-1]
                if pair not in blosum:
                    value = blosum[seq2[j-1],seq1[i-1]]
                else:
                    value = blosum[pair]
                a3 = Wmin[i-1, j-1] + value
                
                Wmin[i,j] = max([a1, a2, a3])
                
                # Параллельно записываем в специальные матрицы сдвигов значение
                # сдвига по горизонтали и вертикали (то есть в каждой следующей
                # ячейке записываем информацию о том, насколько надо изменить 
                # индексы, чтобы придти в предыдущую ячейку).                  
                                                
                if Wmin[i,j] == a1:
                    FromV[i,j] = 0
                    FromG[i,j] = -1
                elif Wmin[i,j] == a2:
                    FromV[i,j] = -1
                    FromG[i,j] = 0
                else:
                    FromV[i,j] = -1
                    FromG[i,j] = -1
                FromG[i,0] = 0
                FromG[0,j] = -1
                FromV[i,0] = -1
                FromV[0,j] = 0
    
    # Восстанавливаем выравнивание с помощью матриц сдвигов. Если сдвиг был
    # диагональным, это либо match, либо mismatch; если горизонтальный/
    # вертикальный - gap, и тогда в соответственной выравненой строке появится 
    # пропуск.
    
    i = n
    j = m         
    while (i == j == 0) == False :
        
        a = FromV[i,j]
        b = FromG[i,j] 
    
        if a == b == -1:
            s1 += seq1[i-1]
            s2 += seq2[j-1]
            i+=-1
            j+=-1
            
        elif a == 0:
            s1 += '-'
            s2 += seq2[j-1]
            i+=0
            j+=-1
           
        elif b == 0:
            s1 += seq1[i-1]
            s2 += '-'
            i+=-1
            j+=0
            
    # Стоимость выравнивания в результате будет равна сумме длин ребер, по 
    # которым пройдет кратчайший путь, и окажется в конце пути - в ячейке 
    # матрицы [n,m].
    # Если нужно кроме выравненных последовательностей и скоров вывести матрицу
    # похожести, раскомментируем:
    # print('Матрица похожести:')
    # print(Wmin)
   
    # Преобразуем список в строку для красоты.
    al_seq1 = ''.join(s1[::-1])
    al_seq2 = ''.join(s2[::-1]) 
    if print_flag:
        print(al_seq1)
        print(al_seq2)
    
    Score = Wmin[n,m]
  
    return Score

# Открываем fasta-файл с несколькими последовательностями.

        
if __name__ == '__main__':
    in_file, stringg, gapp, print_flag, name = parse_inputs()
    
    blosum = MatrixInfo.blosum62
    
    with open(in_file, "r") as handle:
    # Формируем для удобства список из последовательностей.
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        # Заполняем нулями матрицу скоров
        score_matrix = np.zeros((len(records), len(records)))
        # Проходим циклом по последовательностям из файла. Каждую 
        # последовательность выравниваем на саму себя и на другие. Записываем 
        # скор в матрицу скоров.
        for i in range(len(records)):
            for j in range(i,len(records)):
                seq1 = str(records[i].seq)
                seq2 = str(records[j].seq)
                score_matrix[i,j] = levenstein(gapp, print_flag, seq1,seq2)
                score_matrix[j,i] = score_matrix[i,j]
        if print_flag:        
            print(stringg)
            print(score_matrix)
    with open('./{0}.txt'.format(name), 'w') as  out:
        out.write(stringg)
        out.write('\n')
        out.write(str(score_matrix))
    flag = levenstein(gapp, print_flag, seq1, seq2)
