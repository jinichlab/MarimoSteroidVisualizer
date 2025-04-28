import marimo

__generated_with = "0.10.2"
app = marimo.App()


@app.cell
def simple_ui():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    text = mo.ui.text(label="Enter your name:")
    button = mo.ui.run_button(label="Say Hello")
    return button, text


@app.cell
def _(button):
    button
    return


@app.cell
def _(text):
    text
    return


@app.cell
def _(button, text):
    if button.value:
        print(f"Hello, {text.value}!")
    else:
        print("Please type your name and click the button.")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
