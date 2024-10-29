import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def read_matrix_from_file(file_path):
    with open(file_path, 'r') as file:
        # Читаем файл и преобразуем в список списков
        matrix = [list(map(float, line.split())) for line in file]
    return matrix

def plot_heatmap(matrix):
    # Преобразуем в DataFrame для удобства работы с Seaborn
    df = pd.DataFrame(matrix)
    
    # Создаем тепловую карту
    plt.figure(figsize=(8, 6))
    sns.heatmap(df, annot=True, cmap='viridis')
    plt.title('График ошибок по области')
    plt.show()

if __name__ == "__main__":
    # Задайте путь к файлу с матрицей
    file_path = 'matrix.txt'  # Укажите ваш файл с матрицей
    matrix = read_matrix_from_file(file_path)
    plot_heatmap(matrix)
