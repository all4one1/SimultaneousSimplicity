import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog

def read_field_data(filename):
    # Читаем данные из файла
    data = np.loadtxt(filename, skiprows=1)  # Пропускаем заголовок
    x = data[:, 0]
    y = data[:, 1]
    field = data[:, 3]

    # Находим уникальные значения x и y для создания сетки
    unique_x = np.unique(x)
    unique_y = np.unique(y)
    
    # Создаем сетку
    X, Y = np.meshgrid(unique_x, unique_y)
    
    # Преобразуем данные в формат сетки
    Z = np.zeros_like(X)
    for i, (xi, yi, zi) in enumerate(zip(x, y, field)):
        x_idx = np.where(unique_x == xi)[0][0]
        y_idx = np.where(unique_y == yi)[0][0]
        Z[y_idx, x_idx] = zi
    
    return X, Y, Z

def main():
    # Создаем корневое окно tkinter (но скрываем его)
    root = tk.Tk()
    root.withdraw()

    # Открываем диалог выбора файла
    filename = filedialog.askopenfilename(
        title="Выберите файл с данными",
        filetypes=[("Текстовые файлы", "*.txt"), ("Все файлы", "*.*")]
    )

    if not filename:
        print("Файл не выбран")
        return

    try:
        X, Y, Z = read_field_data(filename)
        
        # Создаем график
        plt.figure(figsize=(10, 8))
        
        # Создаем контурный график
        plt.contourf(X, Y, Z, levels=20, cmap='viridis')
        plt.colorbar(label='Значение поля')
        plt.title('Визуализация 2D поля')
        plt.xlabel('X')
        plt.ylabel('Y')
        
        # Устанавливаем равное соотношение сторон
        plt.gca().set_aspect('equal')
        
        plt.show()
        
    except FileNotFoundError:
        print(f"Ошибка: Файл {filename} не найден")
    except Exception as e:
        print(f"Произошла ошибка при обработке файла: {str(e)}")
        print("Убедитесь, что файл имеет правильный формат:")
        print('"x" "y" "field"')
        print("1.0 1.0 0.5")
        print("1.0 2.0 0.7")
        print("...")

if __name__ == "__main__":
    main()
