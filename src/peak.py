class Peak:
    def __init__(self, position: int, height: float, resolution: int) -> None:
        self.positionInReference = position * resolution
        self.position = position
        self.height = height
