import click
import sys


class Logger:
    """
    custom logging for kSpider, shared by click context in the click_context:cli.
    """

    RED = "red"
    BLUE = "blue"
    GREEN = "green"
    YELLOW = "yellow"

    ACTIVE = True

    def __init__(self, active=True):
        self.ACTIVE = active

    def SUCCESS(self, msg):
        if not self.ACTIVE:
            click.secho(f"[SUCCESS] {msg}", fg=Logger.GREEN, bold=True, file=sys.stderr)

    def INFO(self, msg):
        if not self.ACTIVE:
            click.secho(f"[INFO] {msg}", fg=Logger.YELLOW, bold=True, file=sys.stderr)

    def WARNING(self, msg):
        if not self.ACTIVE:
            click.secho(f"[WARNING] {msg}", fg=Logger.YELLOW, bold=True, file=sys.stderr)

    def ERROR(self, msg):        
        click.secho(f"[ERROR] {msg}", fg=Logger.RED, bold=True, file=sys.stderr)
        sys.exit(1)