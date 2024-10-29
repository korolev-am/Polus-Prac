import matplotlib.pyplot as plt

def plot_graph(x_values, y_values):
    # Проверяем, что массивы одинаковой длины
    if len(x_values) != len(y_values):
        raise ValueError("Оба массива должны быть одинаковой длины.")
    
    # Создаем график
    plt.figure(figsize=(10, 5))
    plt.plot(x_values, y_values, marker='o', linestyle='-', color='b')
    
    # Добавляем заголовок и метки осей
    plt.title('График ошибок по итерациям')
    plt.xlabel('Номер итерации')
    plt.ylabel('Величина ошибки')
    
    # Отображаем сетку
    plt.grid()
    
    # Показываем график
    plt.show()

# Пример использования функции
step = 12000
x = [(x + 1)*step for x in range(-1, 17)]
y = [0.000155, 0.000115, 0.000088, 0.000067, 0.000051, 0.000039, 0.000030, 0.000023, 0.000017, 0.000013, 0.000010, 0.000008, 0.000006, 0.000005, 0.000004, 0.000003, 0.000002, 0.000002]
plot_graph(x, y)