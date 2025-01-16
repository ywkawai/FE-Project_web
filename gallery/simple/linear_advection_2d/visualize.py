import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation

u = xr.open_mfdataset("history.pe000000.nc", decode_times=False, combine='by_coords')["u"]
x = u.coords["x"]
y = u.coords["y"]
time = u.coords["time"]

fig, ax = plt.subplots(figsize=(6,4))
pcl = ax.pcolormesh( x, y, u.isel(time=0), vmin=-1.05, vmax=1.05, cmap="jet")
fig.colorbar(pcl, ax=ax)

title = ax.set_title("time=0")

ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("2D linear advection") 

def update(n):
  print(f"n={n}")
  pcl.set_array(u.isel(time=n).values.ravel())
  time_txt = "{:.2f}".format(time.values[n])
  title.set_text("time="+time_txt)

ani = animation.FuncAnimation(fig, update, frames=len(time), interval=200)
ani.save("advect2D.mp4", writer="ffmpeg")