import tkinter as tk
from PIL import Image, ImageTk
import json
import numpy as np

class AdvancedCoordinatePicker:
    def __init__(self, image_path, x_range=(0, 100), y_range=(0, 100)):
        self.root = tk.Tk()
        self.root.title("Advanced Coordinate Picker")
        
        self.x_range = x_range
        self.y_range = y_range

        self.original_image = Image.open(image_path)
        self.image = self.original_image.copy()
        self.photo = ImageTk.PhotoImage(self.image)

        self.scale_factor = 1.0

        self.canvas_frame = tk.Frame(self.root)
        self.canvas_frame.pack(fill=tk.BOTH, expand=True)
        
        self.canvas = tk.Canvas(self.canvas_frame)
        self.scroll_y = tk.Scrollbar(self.canvas_frame, orient="vertical", command=self.canvas.yview)
        self.scroll_x = tk.Scrollbar(self.root, orient="horizontal", command=self.canvas.xview)
        
        self.canvas.configure(yscrollcommand=self.scroll_y.set, 
                             xscrollcommand=self.scroll_x.set)
        
        self.scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        self.scroll_x.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.image_container = self.canvas.create_image(0, 0, anchor=tk.NW, image=self.photo)
        
        self.pixel_coords = []
        self.real_coords = []
        
        self.canvas.bind("<Button-1>", self.on_click)
        self.canvas.bind("<Button-3>", self.on_right_click)
        self.canvas.bind("<MouseWheel>", self.on_mousewheel)
        self.canvas.bind("<Button-4>", self.on_mousewheel)
        self.canvas.bind("<Button-5>", self.on_mousewheel)
        
        self.control_frame = tk.Frame(self.root)
        self.control_frame.pack(fill=tk.X)
        
        tk.Label(self.control_frame, text="X диапазон:").pack(side=tk.LEFT)
        self.x_min_entry = tk.Entry(self.control_frame, width=8)
        self.x_min_entry.insert(0, str(x_range[0]))
        self.x_min_entry.pack(side=tk.LEFT)
        
        tk.Label(self.control_frame, text="до").pack(side=tk.LEFT)
        self.x_max_entry = tk.Entry(self.control_frame, width=8)
        self.x_max_entry.insert(0, str(x_range[1]))
        self.x_max_entry.pack(side=tk.LEFT)
        
        tk.Label(self.control_frame, text="Y диапазон:").pack(side=tk.LEFT)
        self.y_min_entry = tk.Entry(self.control_frame, width=8)
        self.y_min_entry.insert(0, str(y_range[0]))
        self.y_min_entry.pack(side=tk.LEFT)
        
        tk.Label(self.control_frame, text="до").pack(side=tk.LEFT)
        self.y_max_entry = tk.Entry(self.control_frame, width=8)
        self.y_max_entry.insert(0, str(y_range[1]))
        self.y_max_entry.pack(side=tk.LEFT)
        
        tk.Button(self.control_frame, text="Обновить диапазоны", 
                 command=self.update_ranges).pack(side=tk.LEFT, padx=5)
        tk.Button(self.control_frame, text="Сохранить", 
                 command=self.save_coordinates).pack(side=tk.LEFT)
        tk.Button(self.control_frame, text="Очистить", 
                 command=self.clear_coordinates).pack(side=tk.LEFT)
        tk.Button(self.control_frame, text="Зум +", 
                 command=lambda: self.zoom(1.2)).pack(side=tk.LEFT)
        tk.Button(self.control_frame, text="Зум -", 
                 command=lambda: self.zoom(0.8)).pack(side=tk.LEFT)
        
        self.status = tk.Label(self.root, text="Кликайте по изображению для добавления точек")
        self.status.pack()
        
        self.canvas.config(scrollregion=self.canvas.bbox("all"))
        
    def pixel_to_real(self, x, y):
        """Преобразование пиксельных координат в реальные значения"""
        norm_x = x / self.image.width
        norm_y = y / self.image.height
        
        real_x = self.x_range[0] + norm_x * (self.x_range[1] - self.x_range[0])
        real_y = self.y_range[1] - norm_y * (self.y_range[1] - self.y_range[0])
        
        return real_x, real_y
    
    def on_click(self, event):
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        
        real_x, real_y = self.pixel_to_real(x, y)
        
        self.pixel_coords.append((x, y))
        self.real_coords.append((real_x, real_y))
        
        radius = 5
        self.canvas.create_oval(x-radius, y-radius, x+radius, y+radius, 
                               fill="red", outline="red", tags="point")
        
        label = f"({real_x:.2f}, {real_y:.2f})"
        self.canvas.create_text(x, y-15, text=label, 
                               fill="yellow", font=("Arial", 9), tags="label")
        
        self.status.config(text=f"Точка {len(self.real_coords)}: {label}")
        
    def on_right_click(self, event):
        if self.real_coords:
            self.pixel_coords.pop()
            self.real_coords.pop()
            self.redraw_image()
            self.status.config(text=f"Удалена последняя точка. Осталось: {len(self.real_coords)}")
    
    def on_mousewheel(self, event):
        scale = 1.1
        if event.num == 5 or event.delta < 0:
            scale = 0.9
        
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        
        self.zoom(scale, x, y)
    
    def zoom(self, scale, x=None, y=None):
        self.scale_factor *= scale
        
        new_width = int(self.original_image.width * self.scale_factor)
        new_height = int(self.original_image.height * self.scale_factor)
        
        self.image = self.original_image.resize((new_width, new_height), Image.LANCZOS)
        self.photo = ImageTk.PhotoImage(self.image)
        
        self.canvas.itemconfig(self.image_container, image=self.photo)
        self.canvas.config(scrollregion=self.canvas.bbox("all"))
        
    def update_ranges(self):
        try:
            x_min = float(self.x_min_entry.get())
            x_max = float(self.x_max_entry.get())
            y_min = float(self.y_min_entry.get())
            y_max = float(self.y_max_entry.get())
            
            self.x_range = (x_min, x_max)
            self.y_range = (y_min, y_max)
            
            self.real_coords = []
            for px, py in self.pixel_coords:
                rx, ry = self.pixel_to_real(px, py)
                self.real_coords.append((rx, ry))
            
            self.redraw_image()
            self.status.config(text="Диапазоны обновлены")
        except ValueError:
            self.status.config(text="Ошибка: введите корректные числа")
    
    def redraw_image(self):
        self.canvas.delete("point")
        self.canvas.delete("label")
        
        for i, ((px, py), (rx, ry)) in enumerate(zip(self.pixel_coords, self.real_coords), 1):
            radius = 5
            self.canvas.create_oval(px-radius, py-radius, px+radius, py+radius, 
                                   fill="red", outline="red", tags="point")
            
            label = f"({rx:.2f}, {ry:.2f})"
            self.canvas.create_text(px, py-15, text=label, 
                                   fill="yellow", font=("Arial", 9), tags="label")
    
    def save_coordinates(self):
        data = {
            "image_size": (self.original_image.width, self.original_image.height),
            "x_range": self.x_range,
            "y_range": self.y_range,
            "pixel_coordinates": self.pixel_coords,
            "real_coordinates": self.real_coords,
            "scale_factor": self.scale_factor
        }
        
        with open("coordinates_detailed.json", "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, default=lambda x: float(x) if isinstance(x, np.float32) else x)
        
        with open("coordinates.csv", "w", encoding="utf-8") as f:
            f.write("Index,Pixel_X,Pixel_Y,Real_X,Real_Y\n")
            for i, ((px, py), (rx, ry)) in enumerate(zip(self.pixel_coords, self.real_coords), 1):
                f.write(f"{i},{px},{py},{rx:.4f},{ry:.4f}\n")
        
        self.status.config(text=f"Данные сохранены (JSON, CSV). Точек: {len(self.real_coords)}")
        
    def clear_coordinates(self):
        self.pixel_coords = []
        self.real_coords = []
        self.redraw_image()
        self.status.config(text="Все точки удалены")
    
    def run(self):
        self.root.mainloop()

if __name__ == "__main__":
    image_path = "image rnkt.bmp"

    app = AdvancedCoordinatePicker(
        image_path, 
        x_range=(0, 10),
        y_range=(0, 100)
    )
    app.run()